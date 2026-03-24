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


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Interactive user input
# ─────────────────────────────────────────────────────────────────────────────

def get_user_inputs():
    """
    Prompt the user interactively for protein family, taxonomic group,
    and output directory name.

    Validates that protein family and taxon are non-empty strings.
    Defaults output directory to "pipeline_output" if left blank.

    Returns
    -------
    tuple : (protein_family str, taxon str, outdir str)
    """
    print("=" * 62)
    print("  Protein Family Conservation & Motif Pipeline")
    print("=" * 62)
    print()
    print("This programme will:")
    print("  1. Fetch protein sequences from NCBI for your query.")
    print("  2. Report species diversity and ask whether to proceed.")
    print("  3. Align sequences and plot conservation.")
    print("  4. Scan for PROSITE motifs (patmatmotifs).")
    print("  5. Optionally run pepstats (physicochemical properties).")
    print()
    print(f"  Note: up to {MAX_SEQS_FETCH} sequences will be fetched from NCBI.")
    print()

    # Get protein family name - must not be empty
    protein_family = input(
        "Enter protein family name\n"
        "  (e.g. glucose-6-phosphatase, ABC transporter, kinase): "
    ).strip()
    if not protein_family:
        print("[ERROR] Protein family name cannot be empty. Exiting.")
        sys.exit(1)

    # Get taxonomic group - must not be empty
    taxon = input(
        "\nEnter taxonomic group (scientific name)\n"
        "  (e.g. Aves, Mammalia, Rodentia, Vertebrata): "
    ).strip()
    if not taxon:
        print("[ERROR] Taxonomic group cannot be empty. Exiting.")
        sys.exit(1)

    # Get output directory - default if blank
    outdir = input(
        "\nEnter output directory name [default: pipeline_output]: "
    ).strip()
    if not outdir:
        outdir = "pipeline_output"

    print()
    return protein_family, taxon, outdir


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Create output directory
# ─────────────────────────────────────────────────────────────────────────────

def make_outdir(outdir):
    """
    Create the output directory. If it already exists, warn the user
    and ask whether to continue (existing files may be overwritten).

    Exits cleanly if the user declines to continue.

    Parameters
    ----------
    outdir : str - path to the output directory
    """
    if os.path.exists(outdir):
        print(f"[WARNING] Output directory '{outdir}' already exists.")
        print("          Existing files may be overwritten.")
        ans = input("          Continue anyway? [y/N]: ").strip().lower()
        if ans != "y":
            print("Exiting. Choose a different output directory and re-run.")
            sys.exit(0)
    else:
        os.makedirs(outdir)
        print(f"[OK] Created output directory: {outdir}\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Search NCBI and fetch sequences via EDirect
# ─────────────────────────────────────────────────────────────────────────────

def fetch_sequences(protein_family, taxon, outdir):
    """
    Use NCBI EDirect (esearch + efetch) to retrieve protein sequences
    in FASTA format.

    NCBI query constructed as:
        "<protein_family>"[Protein Name] AND "<taxon>"[Organism]

    Runs esearch first to count hits before fetching, so the user
    can be warned if the dataset is very large.

    Saves results to <outdir>/sequences_raw.fasta.

    Parameters
    ----------
    protein_family : str
    taxon          : str
    outdir         : str

    Returns
    -------
    tuple : (raw_fasta_path str, query_string str)
    """
    query = f'"{protein_family}"[Protein Name] AND "{taxon}"[Organism]'
    print(f"[Step 1] Searching NCBI protein database...")
    print(f"         Query: {query}\n")

    # First call: count how many sequences match the query
    search_xml = run_cmd(
        ["esearch", "-db", "protein", "-query", query],
        description="esearch (count)"
    )

    # Parse the <Count> tag from the esearch XML output
    count_match = re.search(r"<Count>(\d+)</Count>", search_xml)
    if not count_match:
        print("[ERROR] Could not parse hit count from esearch XML output.")
        print("        Raw output preview:")
        print(search_xml[:400])
        sys.exit(1)

    count = int(count_match.group(1))
    print(f"         NCBI reports {count} sequence(s) matching your query.")

    # Exit with helpful suggestions if nothing is found
    if count == 0:
        print("\n[ERROR] No sequences found for this query.")
        print("        Suggestions:")
        print("          - Check the spelling of the protein family name.")
        print("          - Try a broader taxonomic group (e.g. Vertebrata).")
        print("          - Try a simpler protein name (e.g. 'kinase').")
        sys.exit(1)

    # Warn if dataset is very large and ask for confirmation
    if count > MAX_SEQS_FETCH:
        print(f"\n[WARNING] {count} sequences found; limit is {MAX_SEQS_FETCH}.")
        print(f"          Only the first {MAX_SEQS_FETCH} will be retrieved.")
        print("          Consider using a more specific taxon or protein name.")
        ans = input(f"          Proceed with {MAX_SEQS_FETCH} sequences? [y/N]: "
                    ).strip().lower()
        if ans != "y":
            print("Exiting.")
            sys.exit(0)
        retmax = MAX_SEQS_FETCH
    else:
        retmax = count

    # Second call: fetch sequences in FASTA format
    # esearch XML output is piped into efetch via input_data argument
    print(f"\n[Step 1] Fetching {retmax} sequence(s) from NCBI...")

    esearch_xml = run_cmd(
        ["esearch", "-db", "protein", "-query", query],
        description="esearch (for fetch)"
    )
    fasta_data = run_cmd(
        ["efetch", "-db", "protein", "-format", "fasta",
         "-retmax", str(retmax)],
        description="efetch",
        input_data=esearch_xml
    )

    # Check efetch actually returned something
    if not fasta_data.strip():
        print("[ERROR] efetch returned no data.")
        print("        Check your internet connection and NCBI availability.")
        sys.exit(1)

    # Save raw sequences to file
    raw_fasta = os.path.join(outdir, "sequences_raw.fasta")
    with open(raw_fasta, "w") as fh:
        fh.write(fasta_data)

    print(f"[OK] Raw sequences saved to: {raw_fasta}")
    return raw_fasta, query
