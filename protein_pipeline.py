#!/usr/bin/env python3
"""
protein_pipeline.py
====================
A generic bioinformatics pipeline for protein family analysis.

Given a user-defined protein family and taxonomic group, this
programme will:
  1. Fetch matching protein sequences from NCBI
  2. Assess species diversity in the retrieved sequences
  3. Perform multiple sequence alignment and conservation analysis
  4. Scan sequences for PROSITE motifs
  5. Optionally run physicochemical analysis with pepstats

Usage:
    python3 protein_pipeline.py

Dependencies (must be on PATH):
    esearch, efetch  (NCBI EDirect)
    clustalo         (Clustal Omega)
    plotcon          (EMBOSS)
    patmatmotifs     (EMBOSS)
    pepstats         (EMBOSS)

Author: Rahul Dey
"""

import subprocess
import sys
import os
import re
import collections


# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

# Maximum number of sequences to use for alignment and conservation plotting.
# If more sequences than this are retrieved, a representative subset is chosen.
MAX_SEQS_CONSERVATION = 100

# If NCBI returns more hits than this, warn the user and ask for confirmation
# before proceeding. Not a hard block - user can choose to continue.
MAX_SEQS_FETCH = 1000


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: safely run an external command and return stdout
# ─────────────────────────────────────────────────────────────────────────────

def run_cmd(cmd, description="command", input_data=None):
    """
    Run a shell command given as a list of strings.
    Returns stdout as a string.
    Exits with an informative message on failure.

    All external tool calls in this pipeline go through this function
    so that errors are handled consistently in one place.

    Parameters
    ----------
    cmd         : list of str  - command + arguments
    description : str          - human-readable label for error messages
    input_data  : str or None  - optional text to pipe into stdin
                                 (used for esearch | efetch chaining)

    Returns
    -------
    str - stdout from the command
    """
    print(f"  [Running] {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            input=input_data,
            capture_output=True,
            text=True
        )
    except FileNotFoundError:
        # Tool not found on PATH at all
        print(f"\n[ERROR] Program '{cmd[0]}' not found on PATH.")
        print("        Please ensure it is installed and accessible.")
        sys.exit(1)

    # Non-zero exit code means the tool reported an error
    if result.returncode != 0:
        print(f"\n[ERROR] {description} failed (exit code {result.returncode}).")
        print(f"        STDERR: {result.stderr.strip()}")
        sys.exit(1)

    return result.stdout


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: verify all required external tools are available before starting
# ─────────────────────────────────────────────────────────────────────────────

def check_dependencies(tools):
    """
    Check each tool name in 'tools' is reachable via the `which` command.

    Collects all missing tools into a list before reporting, so the user
    sees everything that is missing in one go rather than one at a time.

    Exits with a clear message if any tools are absent.

    Parameters
    ----------
    tools : list of str - tool names to check (e.g. ["clustalo", "plotcon"])
    """
    missing = []

    for tool in tools:
        # `which` returns exit code 0 if the tool is found, 1 if not
        result = subprocess.run(
            ["which", tool],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            missing.append(tool)

    if missing:
        print("\n[ERROR] The following required tools were not found on PATH:")
        for m in missing:
            print(f"        - {m}")
        print("        Please install them or add them to PATH before running.")
        sys.exit(1)

    print("[OK] All required tools are available.\n")
