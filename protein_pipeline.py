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
        print(f"\n[ERROR] Program '{cmd[0]}' not found on PATH.")
        print("        Please ensure it is installed and accessible.")
        sys.exit(1)

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

    protein_family = input(
        "Enter protein family name\n"
        "  (e.g. glucose-6-phosphatase, ABC transporter, kinase): "
    ).strip()
    if not protein_family:
        print("[ERROR] Protein family name cannot be empty. Exiting.")
        sys.exit(1)

    taxon = input(
        "\nEnter taxonomic group (scientific name)\n"
        "  (e.g. Aves, Mammalia, Rodentia, Vertebrata): "
    ).strip()
    if not taxon:
        print("[ERROR] Taxonomic group cannot be empty. Exiting.")
        sys.exit(1)

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

    search_xml = run_cmd(
        ["esearch", "-db", "protein", "-query", query],
        description="esearch (count)"
    )

    count_match = re.search(r"<Count>(\d+)</Count>", search_xml)
    if not count_match:
        print("[ERROR] Could not parse hit count from esearch XML output.")
        print("        Raw output preview:")
        print(search_xml[:400])
        sys.exit(1)

    count = int(count_match.group(1))
    print(f"         NCBI reports {count} sequence(s) matching your query.")

    if count == 0:
        print("\n[ERROR] No sequences found for this query.")
        print("        Suggestions:")
        print("          - Check the spelling of the protein family name.")
        print("          - Try a broader taxonomic group (e.g. Vertebrata).")
        print("          - Try a simpler protein name (e.g. 'kinase').")
        sys.exit(1)

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

    if not fasta_data.strip():
        print("[ERROR] efetch returned no data.")
        print("        Check your internet connection and NCBI availability.")
        sys.exit(1)

    raw_fasta = os.path.join(outdir, "sequences_raw.fasta")
    with open(raw_fasta, "w") as fh:
        fh.write(fasta_data)

    print(f"[OK] Raw sequences saved to: {raw_fasta}")
    return raw_fasta, query


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: pure-Python FASTA parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_fasta(filepath):
    """
    Parse a FASTA file into a list of (header, sequence) tuples.
    Header string does NOT include the leading '>'.
    Blank lines within records are silently ignored.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    list of (str, str)
    """
    records   = []
    header    = None
    seq_parts = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header    = line[1:]
                seq_parts = []
            elif line:
                seq_parts.append(line.upper())

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def write_fasta(records, filepath):
    """
    Write a list of (header, sequence) tuples to a FASTA file.
    Sequences are line-wrapped at 60 characters (standard FASTA format).

    Parameters
    ----------
    records  : list of (str, str)
    filepath : str
    """
    with open(filepath, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: extract species name from an NCBI FASTA header
# ─────────────────────────────────────────────────────────────────────────────

def extract_species(header):
    """
    Extract organism name from a standard NCBI FASTA header.
    NCBI headers typically end with [Organism Name] in square brackets.

    Example:
      "XP_012345.1 glucose-6-phosphatase [Gallus gallus]"
      returns "Gallus gallus"

    Falls back to the first two whitespace-delimited tokens if no
    square brackets are present.

    Parameters
    ----------
    header : str

    Returns
    -------
    str - species/organism name
    """
    match = re.search(r"\[([^\]]+)\]", header)
    if match:
        return match.group(1)
    words = header.split()
    return " ".join(words[:2]) if len(words) >= 2 else header


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Species diversity assessment
# ─────────────────────────────────────────────────────────────────────────────

def assess_diversity(records, outdir):
    """
    Tally sequences per species using collections.Counter.
    Print a summary table to screen and save the full table to file.
    Warn if diversity is low. Ask the user whether to continue.

    Parameters
    ----------
    records : list of (str, str) - parsed FASTA records
    outdir  : str

    Returns
    -------
    list of (str, str) - same records passed through unchanged
    """
    print("[Step 2] Assessing species diversity...\n")

    species_counts = collections.Counter(
        extract_species(header) for header, _ in records
    )

    n_seqs    = len(records)
    n_species = len(species_counts)

    print(f"         Total sequences : {n_seqs}")
    print(f"         Unique species  : {n_species}\n")
    print(f"         {'Species':<45} {'Count':>6}")
    print("         " + "-" * 53)
    for sp, cnt in species_counts.most_common(20):
        print(f"         {sp:<45} {cnt:>6}")
    if n_species > 20:
        print(f"         ... and {n_species - 20} more species (see file).")
    print()

    div_file = os.path.join(outdir, "species_diversity.txt")
    with open(div_file, "w") as fh:
        fh.write(f"Total sequences : {n_seqs}\n")
        fh.write(f"Unique species  : {n_species}\n\n")
        fh.write(f"{'Species':<50} {'Count':>6}\n")
        fh.write("-" * 58 + "\n")
        for sp, cnt in species_counts.most_common():
            fh.write(f"{sp:<50} {cnt:>6}\n")
    print(f"[OK] Full diversity report saved to: {div_file}")

    if n_species == 1:
        print("\n[WARNING] All sequences are from a single species.")
        print("          Conservation analysis across species will not be meaningful.")
    elif n_species < 5:
        print(f"\n[WARNING] Only {n_species} species are represented.")
        print("          Consider a broader taxonomic group for richer results.")

    ans = input("\nContinue with this dataset? [y/N]: ").strip().lower()
    if ans != "y":
        print("Exiting. Try a different protein family or taxonomic group.")
        sys.exit(0)

    return records


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Subsample sequences for alignment if needed
# ─────────────────────────────────────────────────────────────────────────────

def subsample_for_alignment(records, outdir):
    """
    If the number of records exceeds MAX_SEQS_CONSERVATION, select a
    representative subset by:
      1. Sorting sequences by descending length (favouring complete sequences)
      2. Picking evenly-spaced indices across the sorted list

    Saves the result to <outdir>/sequences_for_alignment.fasta

    Parameters
    ----------
    records : list of (str, str)
    outdir  : str

    Returns
    -------
    tuple : (align_fasta_path str, n_seqs int)
    """
    n = len(records)
    align_fasta = os.path.join(outdir, "sequences_for_alignment.fasta")

    if n <= MAX_SEQS_CONSERVATION:
        write_fasta(records, align_fasta)
        print(f"\n[Step 3] Using all {n} sequences for alignment.")
    else:
        print(f"\n[Step 3] {n} sequences found; subsampling to"
              f" {MAX_SEQS_CONSERVATION} for alignment.")
        print("         Method: sort by descending sequence length,"
              " pick evenly-spaced subset.")

        sorted_recs = sorted(records, key=lambda r: len(r[1]), reverse=True)
        step        = n / MAX_SEQS_CONSERVATION
        indices     = [int(i * step) for i in range(MAX_SEQS_CONSERVATION)]
        subset      = [sorted_recs[i] for i in indices]

        write_fasta(subset, align_fasta)
        print(f"[OK] Subsampled FASTA ({MAX_SEQS_CONSERVATION} seqs) saved to:"
              f" {align_fasta}")

    return align_fasta, min(n, MAX_SEQS_CONSERVATION)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Multiple sequence alignment with clustalo
# ─────────────────────────────────────────────────────────────────────────────

def run_alignment(input_fasta, outdir):
    """
    Align sequences using Clustal Omega (clustalo).
    Output format is Pearson/FASTA (--outfmt=fasta).

    Parameters
    ----------
    input_fasta : str - path to unaligned FASTA
    outdir      : str

    Returns
    -------
    str - path to the aligned FASTA file
    """
    aligned_fasta = os.path.join(outdir, "sequences_aligned.fasta")
    print(f"\n[Step 4] Running multiple sequence alignment (clustalo)...")

    run_cmd(
        [
            "clustalo",
            "-i",        input_fasta,
            "-o",        aligned_fasta,
            "--outfmt=fasta",
            "--force",
            "-v"
        ],
        description="clustalo"
    )

    print(f"[OK] Aligned FASTA saved to: {aligned_fasta}")
    return aligned_fasta


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Conservation analysis - EMBOSS plotcon + Python column scoring
# ─────────────────────────────────────────────────────────────────────────────

def run_plotcon(aligned_fasta, outdir):
    """
    (a) Run EMBOSS plotcon to produce a graphical conservation plot (PNG).
    (b) Compute a pure-Python per-column identity score from the alignment,
        save as a text file, and print an ASCII profile to screen.

    Parameters
    ----------
    aligned_fasta : str - path to aligned FASTA
    outdir        : str
    """
    print(f"\n[Step 5] Running EMBOSS plotcon (conservation plot)...")

    goutfile_base = os.path.join(outdir, "conservation_plot")

    run_cmd(
        [
            "plotcon",
            "-sequences", aligned_fasta,
            "-winsize",   "4",
            "-graph",     "png",
            "-goutfile",  goutfile_base
        ],
        description="plotcon"
    )

    # EMBOSS plotcon appends .1.png to the goutfile name
    emboss_png = goutfile_base + ".1.png"
    target_png = goutfile_base + ".png"

    if os.path.exists(emboss_png):
        os.rename(emboss_png, target_png)
        print(f"[OK] Conservation plot saved to: {target_png}")
    elif os.path.exists(target_png):
        print(f"[OK] Conservation plot saved to: {target_png}")
    else:
        print("[WARNING] plotcon PNG not found at expected path.")
        print("          Check the output directory for any .png files.")

    # Run the pure-Python conservation scoring as well
    compute_python_conservation(aligned_fasta, outdir)


def compute_python_conservation(aligned_fasta, outdir):
    """
    Pure-Python per-column conservation scoring for a multiple alignment.

    For each alignment column, calculates the fraction of non-gap residues
    that match the most common residue (identity-based conservation score).
    Score range: 0.0 (fully variable) to 1.0 (fully identical).

    Saves per-column scores to a tab-delimited text file.
    Prints a summary and ASCII bar chart to screen.

    Parameters
    ----------
    aligned_fasta : str - path to aligned FASTA
    outdir        : str
    """
    records = parse_fasta(aligned_fasta)
    if not records:
        print("[WARNING] Could not parse aligned FASTA; skipping conservation.")
        return

    sequences = [seq.upper() for _, seq in records]
    n_seqs    = len(sequences)
    aln_len   = max(len(s) for s in sequences)

    # Pad any shorter sequences to alignment length with gap characters
    sequences = [s.ljust(aln_len, "-") for s in sequences]

    # Score each column using Counter to find the most common residue
    scores = []
    for col in range(aln_len):
        column  = [sequences[i][col] for i in range(n_seqs)]
        non_gap = [c for c in column if c != "-"]
        if not non_gap:
            scores.append(0.0)
            continue
        top_count = collections.Counter(non_gap).most_common(1)[0][1]
        scores.append(top_count / len(non_gap))

    avg_score = sum(scores) / len(scores)
    high_cons = sum(1 for s in scores if s >= 0.8)

    print(f"\n         [Python Conservation Summary]")
    print(f"         Alignment length        : {aln_len} columns")
    print(f"         Number of sequences     : {n_seqs}")
    print(f"         Mean conservation score : {avg_score:.3f}")
    print(f"         Columns >= 80% identical: {high_cons}"
          f" ({100 * high_cons / aln_len:.1f}%)")

    # Save tab-delimited scores file
    cons_file = os.path.join(outdir, "conservation_scores.txt")
    with open(cons_file, "w") as fh:
        fh.write("Column\tConservation_Score\n")
        for i, s in enumerate(scores, 1):
            fh.write(f"{i}\t{s:.4f}\n")
    print(f"[OK] Per-column conservation scores saved to: {cons_file}")

    # Print ASCII bar chart sampled across the alignment
    print("\n         Conservation profile (sampled every ~2%):")
    sample_step = max(1, aln_len // 50)
    for col in range(0, aln_len, sample_step):
        bar = "█" * int(scores[col] * 30)
        print(f"         {col+1:>6} | {bar:<30}  {scores[col]:.2f}")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: PROSITE motif scan with EMBOSS patmatmotifs
# ─────────────────────────────────────────────────────────────────────────────

def parse_patmatmotifs(result_file):
    """
    Parse a patmatmotifs output file and return a list of unique motif labels.

    patmatmotifs output contains lines like:
        Motif = PS00107  PROTEIN_KINASE_ATP

    Parameters
    ----------
    result_file : str - path to patmatmotifs output file

    Returns
    -------
    list of str - unique motif labels found
    """
    motifs = set()
    if not os.path.exists(result_file):
        return []

    with open(result_file, "r") as fh:
        for line in fh:
            m = re.search(r"Motif\s*=\s*(\S+)\s*(.*)", line)
            if m:
                motif_id   = m.group(1).strip()
                motif_name = m.group(2).strip()
                motifs.add(f"{motif_id} {motif_name}".strip())

    return list(motifs)


def run_prosite_scan(records, outdir):
    """
    Run EMBOSS patmatmotifs on each sequence individually against the
    PROSITE motif database.

    Individual results saved under <outdir>/prosite_results/.
    Combined summary saved to <outdir>/prosite_summary.txt.

    Note: patmatmotifs is called via subprocess.run() directly rather
    than run_cmd() because some EMBOSS versions return a non-zero exit
    code when no motifs are found, which is not an error condition.

    Parameters
    ----------
    records : list of (str, str)
    outdir  : str
    """
    print(f"\n[Step 6] Scanning for PROSITE motifs (patmatmotifs)...")

    prosite_dir = os.path.join(outdir, "prosite_results")
    os.makedirs(prosite_dir, exist_ok=True)

    # motif_label -> list of accession strings that contain the motif
    motif_summary = collections.defaultdict(list)
    n_with_hits   = 0

    for i, (header, seq) in enumerate(records):

        # Derive a safe filename from the first token of the header
        acc      = header.split()[0].replace("/", "_").replace("|", "_")
        seq_file = os.path.join(prosite_dir, f"{acc}.fasta")
        out_file = os.path.join(prosite_dir, f"{acc}.patmatmotifs")

        # Write single-sequence FASTA for patmatmotifs input
        with open(seq_file, "w") as fh:
            fh.write(f">{header}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j+60] + "\n")

        # Run patmatmotifs - non-fatal if non-zero exit (no hits case)
        subprocess.run(
            ["patmatmotifs", "-sequence", seq_file,
             "-outfile", out_file, "-full", "Y"],
            capture_output=True,
            text=True
        )

        # Parse any motifs found from the output file
        motifs = parse_patmatmotifs(out_file)
        if motifs:
            n_with_hits += 1
            for motif in motifs:
                motif_summary[motif].append(acc)

        # Print a progress counter every 10 sequences
        if (i + 1) % 10 == 0 or (i + 1) == len(records):
            print(f"         Scanned {i+1}/{len(records)} sequences...", end="\r")

    print()  # newline after progress line

    # Write the combined motif summary report
    summary_file = os.path.join(outdir, "prosite_summary.txt")
    with open(summary_file, "w") as fh:
        fh.write("PROSITE Motif Scan Summary\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Sequences scanned      : {len(records)}\n")
        fh.write(f"Sequences with hits    : {n_with_hits}\n\n")
        if motif_summary:
            fh.write(f"{'Motif':<45} {'Occurrences':>12}\n")
            fh.write("-" * 59 + "\n")
            for motif, accs in sorted(motif_summary.items(),
                                      key=lambda x: -len(x[1])):
                fh.write(f"{motif:<45} {len(accs):>12}\n")
        else:
            fh.write("No PROSITE motifs found in any sequence.\n")

    print(f"\n[OK] PROSITE scan complete.")
    print(f"     Sequences with hits : {n_with_hits}/{len(records)}")
    print(f"     Unique motifs found : {len(motif_summary)}")
    if motif_summary:
        print("\n     Top motifs:")
        for motif, accs in sorted(motif_summary.items(),
                                  key=lambda x: -len(x[1]))[:10]:
            print(f"       {motif:<45} in {len(accs)} sequence(s)")
    print(f"[OK] Summary saved to: {summary_file}")
