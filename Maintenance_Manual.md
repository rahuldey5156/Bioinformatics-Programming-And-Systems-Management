================================================================
  Protein Family Conservation & Motif Pipeline
  USER MANUAL — Section 2: Maintenance Manual for Programmers
================================================================

OVERVIEW
---------
protein_pipeline.py is a single-file Python 3 script (~500 lines)
implementing a linear bioinformatics pipeline. It uses only the
Python standard library (subprocess, sys, os, re, collections).
No BioPython. No third-party packages.

External tools are called via subprocess.run() through a single
wrapper function (run_cmd). All file I/O uses built-in open().
FASTA parsing is done in pure Python.

The script is structured as a set of discrete functions, each
corresponding to one pipeline step, all orchestrated by main().

----------------------------------------------------------------
DEPENDENCIES
----------------------------------------------------------------
Python >= 3.6  (uses f-strings and subprocess capture_output)

External tools (must be on PATH):
  esearch        NCBI EDirect - database search
  efetch         NCBI EDirect - sequence retrieval
  clustalo       Clustal Omega - multiple sequence alignment
  plotcon        EMBOSS - conservation plot generation
  patmatmotifs   EMBOSS - PROSITE motif scanning
  pepstats       EMBOSS - physicochemical properties

----------------------------------------------------------------
CONSTANTS (top of file)
----------------------------------------------------------------
MAX_SEQS_CONSERVATION = 100
  Hard cap on sequences passed to clustalo and plotcon.
  Above this, subsample_for_alignment() reduces the set.
  Increase if the server can handle longer alignment runtimes.

MAX_SEQS_FETCH = 1000
  Threshold at which the user is warned and asked to confirm
  before fetching from NCBI. Not a hard block - user can
  proceed. Increase or decrease as appropriate.

----------------------------------------------------------------
FUNCTION REFERENCE
----------------------------------------------------------------

run_cmd(cmd, description, input_data=None)
------------------------------------------
Central wrapper for all external subprocess calls.

  Arguments:
    cmd         : list of str  - command + arguments
    description : str          - label for error messages
    input_data  : str or None  - piped stdin (used for
                                 esearch | efetch chaining)
  Returns:
    str - stdout of the command

  Behaviour:
    - Calls subprocess.run() with capture_output=True, text=True
    - Catches FileNotFoundError -> prints clean message, exits
    - Checks returncode != 0   -> prints stderr, exits
    - Returns stdout as string

  Note: patmatmotifs legitimately returns non-zero when it finds
  no hits in some EMBOSS versions. For this reason, the PROSITE
  scan step calls subprocess.run() directly rather than run_cmd(),
  and treats non-zero exit as non-fatal.

check_dependencies(tools)
--------------------------
Called once at the start of main() before any user input.
Iterates over a list of tool name strings, runs `which <tool>`
for each, and collects any that return non-zero into a `missing`
list. If missing is non-empty, prints the full list and exits.
This prevents the pipeline failing halfway through due to a
missing tool.

get_user_inputs()
------------------
Prints a welcome banner and pipeline description, then collects
three strings via input():
  - protein_family  : validated non-empty
  - taxon           : validated non-empty
  - outdir          : defaults to "pipeline_output" if blank

Returns: (protein_family, taxon, outdir) as a tuple of strings.

make_outdir(outdir)
--------------------
Creates the output directory with os.makedirs().
If the directory already exists, warns the user and asks for
confirmation before continuing. Exits cleanly if user says no.

fetch_sequences(protein_family, taxon, outdir)
-----------------------------------------------
Constructs the NCBI query string:
  "<protein_family>"[Protein Name] AND "<taxon>"[Organism]

Runs esearch twice:
  1. First call: gets XML output, parses <Count> tag with regex
     to determine how many hits exist before fetching.
  2. Second call: pipes esearch XML into efetch via input_data
     to retrieve FASTA sequences.

The two-step approach (count then fetch) allows the programme
to warn and seek confirmation before retrieving large datasets.

efetch is called with -retmax set to min(count, MAX_SEQS_FETCH).

Output written to: <outdir>/sequences_raw.fasta
Returns: (raw_fasta path, query string)

parse_fasta(filepath)
----------------------
Pure-Python FASTA parser. Iterates over file lines:
  - Lines starting with ">" set a new header (strip the ">")
  - Other non-blank lines are accumulated into seq_parts list
  - On each new ">" line, the previous (header, seq) is appended
  - Final record appended after the loop

Returns: list of (header str, sequence str) tuples.
Sequences are uppercased at point of use, not here.

write_fasta(records, filepath)
-------------------------------
Writes a list of (header, seq) tuples to a FASTA file.
Sequences are wrapped at 60 characters per line using a
range(0, len(seq), 60) slice loop. Standard FASTA format.

extract_species(header)
------------------------
Extracts organism name from an NCBI FASTA header string.
Uses re.search(r"\[([^\]]+)\]", header) to find text inside
the last square brackets (standard NCBI format).
Falls back to " ".join(header.split()[:2]) if no brackets.

assess_diversity(records, outdir)
----------------------------------
Builds a collections.Counter from extract_species() applied
to every header in records. Prints a formatted table of the
top 20 species to screen. Writes the full table to
species_diversity.txt.

Issues warnings (non-fatal) if:
  n_species == 1  (single species - conservation not meaningful)
  n_species < 5   (low diversity - results may be limited)

Asks user to confirm continuation. Exits cleanly if user
declines. Returns records unchanged (pass-through).

subsample_for_alignment(records, outdir)
-----------------------------------------
If len(records) <= MAX_SEQS_CONSERVATION: writes all records
to sequences_for_alignment.fasta and returns.

Otherwise:
  1. Sorts records by descending sequence length using
     sorted(..., key=lambda r: len(r[1]), reverse=True)
     Rationale: longer sequences are more likely to be complete,
     full-length entries rather than fragments.
  2. Calculates step = n / MAX_SEQS_CONSERVATION
  3. Picks indices [int(i * step) for i in range(MAX_SEQS_CONSERVATION)]
     This gives an evenly-spaced sample across the length-sorted list,
     capturing both long and shorter sequences representatively.

Writes subset to sequences_for_alignment.fasta.
Returns: (align_fasta path, n_seqs int)

run_alignment(input_fasta, outdir)
-----------------------------------
Calls clustalo with:
  -i        input FASTA
  -o        output FASTA (sequences_aligned.fasta)
  --outfmt=fasta
  --force   overwrite output if it exists
  -v        verbose (progress to stderr, captured by run_cmd)

Returns: path to aligned FASTA.

run_plotcon(aligned_fasta, outdir)
-----------------------------------
Calls plotcon with:
  -sequences  aligned FASTA
  -winsize 4  sliding window of 4 residues
  -graph png  output as PNG image
  -goutfile   base path for output file

EMBOSS plotcon appends ".1.png" to the goutfile name.
The function renames this to conservation_plot.png.
If neither path exists, a warning is printed (non-fatal).

Then calls compute_python_conservation() for the Python-based
scoring (see below).

compute_python_conservation(aligned_fasta, outdir)
---------------------------------------------------
Pure-Python per-column conservation scoring.

Algorithm:
  For each column index col in range(alignment_length):
    - Extract column = [seq[col] for seq in sequences]
    - Filter gaps: non_gap = [c for c in column if c != "-"]
    - If non_gap is empty: score = 0.0
    - Else: score = most_common_count / len(non_gap)
      where most_common_count comes from
      collections.Counter(non_gap).most_common(1)[0][1]

This gives a simple identity-based conservation score per
column: 1.0 = all non-gap residues identical, 0.0 = maximally
variable.

Summary statistics computed in Python:
  avg_score = sum(scores) / len(scores)
  high_cons = count of columns where score >= 0.8

Output:
  conservation_scores.txt  (tab-delimited: Column, Score)
  ASCII bar chart printed to screen (sampled every ~2%)

parse_patmatmotifs(result_file)
--------------------------------
Parses a single patmatmotifs output file.
Uses re.search(r"Motif\s*=\s*(\S+)\s*(.*)", line) to find
lines reporting a motif hit. Builds a set of
"PSXXXXX MOTIF_NAME" label strings (deduped per sequence).
Returns list of unique motif label strings.
Returns empty list if the file does not exist.

run_prosite_scan(records, outdir)
----------------------------------
Iterates over all records (not just the alignment subset -
the full fetched set is scanned).

For each record:
  1. Writes a single-sequence FASTA to prosite_results/<acc>.fasta
  2. Calls patmatmotifs via subprocess.run() directly (not
     run_cmd) because some EMBOSS versions return non-zero
     exit when no motifs are found, which is not an error.
  3. Calls parse_patmatmotifs() on the output file.
  4. Appends hits to motif_summary defaultdict(list):
       motif_label -> [acc1, acc2, ...]

Prints a progress counter every 10 sequences (using end="\r"
to overwrite the same line).

Writes prosite_summary.txt: sorted by descending occurrence
count, showing motif label and number of sequences it hit.
Prints top 10 motifs to screen.

run_pepstats(raw_fasta, outdir)
--------------------------------
Calls pepstats on the full sequences_raw.fasta in one call
(pepstats accepts multi-sequence FASTA input).

Parses the output file with two regex patterns:
  re.search(r"Molecular weight\s*=\s*([\d.]+)", line)
  re.search(r"Isoelectric Point\s*=\s*([\d.]+)", line)

Collects all MW and pI values into lists, computes mean/min/max
using pure Python (sum/min/max builtins), and prints the summary.

write_summary_report(...)
--------------------------
Writes pipeline_summary.txt containing:
  - protein_family, taxon, query string
  - n_fetched, n_aligned counts
  - Directory listing of all output files with sizes
    (os.path.getsize for each file in os.listdir(outdir))

main()
-------
Orchestrates all steps in order. No logic of its own beyond
calling functions and passing return values between them.
The only branching in main() is the optional pepstats prompt.

----------------------------------------------------------------
PIPELINE FLOW DIAGRAM
----------------------------------------------------------------

  START
    |
    v
  check_dependencies()
  [exits if any tool missing]
    |
    v
  get_user_inputs()
  --> protein_family, taxon, outdir
    |
    v
  make_outdir()
  [warns + confirms if exists]
    |
    v
  fetch_sequences()
    |-- esearch (count XML) --> regex parse <Count>
    |   [exits if 0 hits]
    |   [warns + confirms if > MAX_SEQS_FETCH]
    |-- esearch | efetch --> sequences_raw.fasta
    --> returns (raw_fasta, query)
    |
    v
  parse_fasta(raw_fasta)
  --> records : list of (header, seq)
    |
    v
  assess_diversity(records)
    |-- Counter of species
    |-- writes species_diversity.txt
    |-- warns if n_species < 5
    |-- user: continue? [y/N]
    |   [exits if no]
    --> returns records (unchanged)
    |
    v
  subsample_for_alignment(records)
    |-- if n > MAX_SEQS_CONSERVATION:
    |     sort by len DESC
    |     pick evenly-spaced indices
    |-- writes sequences_for_alignment.fasta
    --> returns (align_fasta, n_aligned)
    |
    v
  run_alignment(align_fasta)
    |-- clustalo
    |-- writes sequences_aligned.fasta
    --> returns aligned_fasta
    |
    v
  run_plotcon(aligned_fasta)
    |-- plotcon --> conservation_plot.png
    |-- compute_python_conservation()
    |     column-by-column Counter scoring
    |     writes conservation_scores.txt
    |     prints ASCII bar chart
    |
    v
  run_prosite_scan(records)
    |-- for each record:
    |     write temp FASTA
    |     patmatmotifs (subprocess direct)
    |     parse_patmatmotifs() --> motif labels
    |-- writes prosite_results/ (per-sequence files)
    |-- writes prosite_summary.txt
    |
    v
  [optional] run_pepstats(raw_fasta)
    |-- pepstats on full sequence set
    |-- regex parse MW + pI
    |-- writes pepstats_results.txt
    |
    v
  write_summary_report()
    |-- writes pipeline_summary.txt
    |
    v
  END

----------------------------------------------------------------
DATA FLOW BETWEEN FUNCTIONS
----------------------------------------------------------------

  Function                 Key input              Key output
  ─────────────────────────────────────────────────────────────
  fetch_sequences()        protein_family, taxon  raw_fasta path
                                                  query string
  parse_fasta()            raw_fasta path         records list
  assess_diversity()       records list           records list
                                                  (pass-through)
  subsample_for_alignment  records list           align_fasta path
                                                  n_aligned int
  run_alignment()          align_fasta path       aligned_fasta path
  run_plotcon()            aligned_fasta path     (files only)
  run_prosite_scan()       records list           (files only)
  run_pepstats()           raw_fasta path         (files only)
  write_summary_report()   all of the above       (file only)

The central data structure passed between most functions is:

  records : list of (str, str)
    A list of (header, sequence) tuples.
    - header : NCBI FASTA header with ">" stripped
    - sequence : full amino acid sequence, uppercase

  motif_summary : collections.defaultdict(list)
    Local to run_prosite_scan().
    Maps motif label string -> list of accession strings.
    Used to rank and count motif occurrences across all sequences.

  species_counts : collections.Counter
    Local to assess_diversity().
    Maps species name string -> integer count.

----------------------------------------------------------------
ERROR TRAPPING SUMMARY
----------------------------------------------------------------

  Location                   What is trapped
  ─────────────────────────────────────────────────────────────
  check_dependencies()       Missing tools on PATH -> exit
  get_user_inputs()          Empty protein_family -> exit
                             Empty taxon -> exit
  make_outdir()              Directory exists -> warn + confirm
  fetch_sequences()          Zero NCBI hits -> exit with tips
                             Unparseable <Count> XML -> exit
                             Empty efetch response -> exit
                             Count > MAX_SEQS_FETCH -> confirm
  parse_fasta()              (robust: blank lines ignored)
  assess_diversity()         n_species == 1 -> warn (non-fatal)
                             n_species < 5  -> warn (non-fatal)
                             User declines  -> clean exit
  run_alignment()            clustalo non-zero exit -> exit
  run_plotcon()              PNG not at expected path -> warn
  compute_python_conservation  Empty aligned FASTA -> warn+return
  run_prosite_scan()         patmatmotifs non-zero -> non-fatal
                             Output file missing -> empty list
  run_pepstats()             Unparseable output -> warn (non-fatal)
  run_cmd() (all steps)      FileNotFoundError -> exit + message
                             returncode != 0   -> exit + stderr

----------------------------------------------------------------
HOW TO EXTEND THE CODE
----------------------------------------------------------------

Adding a new analysis step:
  1. Write a new function following this pattern:

       def run_myanalysis(records, outdir):
           '''Docstring explaining inputs and outputs.'''
           out_file = os.path.join(outdir, "myanalysis.txt")
           run_cmd(["mytool", "-i", ..., "-o", out_file],
                   description="mytool")
           # parse output in Python, print summary
           # write results to out_file

  2. Add the new tool to the required_tools list in main():

       required_tools = [..., "mytool"]

  3. Call run_myanalysis() in main() at the appropriate point
     in the pipeline sequence.

Changing the subsampling strategy:
  Edit subsample_for_alignment(). The current approach sorts
  by descending length then evenly spaces.
  Alternative strategies to consider:
    - Cluster with cd-hit first, pick one representative
      per cluster (requires cd-hit on PATH)
    - Randomly sample (use random.sample)
    - Keep one sequence per species (requires extract_species)

Changing the NCBI query format:
  Edit the query string in fetch_sequences():
    query = f'"{protein_family}"[Protein Name] AND "{taxon}"[Organism]'
  NCBI field tags that may be useful:
    [Protein Name]  - searches the protein name field
    [Organism]      - restricts to a taxonomic group
    [Gene Name]     - searches by gene symbol
    [All Fields]    - searches all text fields (broader)
  Reference: https://www.ncbi.nlm.nih.gov/books/NBK49540/

----------------------------------------------------------------
TESTING RECOMMENDATIONS
----------------------------------------------------------------
Recommended test set:
  protein_family : glucose-6-phosphatase
  taxon          : Aves
  outdir         : test_g6pase_birds

This is small (~20-50 sequences), covers multiple species,
has a known PROSITE hit (PS00174), and runs in under 2 minutes
on the MSc server.

To test individual functions in isolation without running the
full pipeline, import them at the Python prompt:

  python3
  >>> from protein_pipeline import parse_fasta, extract_species
  >>> records = parse_fasta("test.fasta")
  >>> print(len(records), "records parsed")
  >>> print(extract_species(records[0][0]))

To test the FASTA writer round-trip:
  >>> from protein_pipeline import parse_fasta, write_fasta
  >>> recs = parse_fasta("sequences_raw.fasta")
  >>> write_fasta(recs, "round_trip_test.fasta")
  >>> recs2 = parse_fasta("round_trip_test.fasta")
  >>> print(len(recs) == len(recs2))   # should be True

----------------------------------------------------------------
KNOWN LIMITATIONS
----------------------------------------------------------------
1. The NCBI query uses [Protein Name] which matches exact name
   tokens. Synonyms (e.g. "G6Pase" vs "glucose-6-phosphatase")
   will not be matched unless the user tries multiple terms.

2. patmatmotifs requires the PROSITE database files to be
   installed on the server (typically pointed to by the
   EMBOSS_DATA environment variable). If motif scanning returns
   no results unexpectedly, check that PROSITE data files exist
   at the expected EMBOSS data path.

3. plotcon output filename behaviour (adding ".1.png") may vary
   between EMBOSS versions. The code handles both ".1.png" and
   ".png" but other suffixes would require editing run_plotcon().

4. For very large datasets (>500 sequences), the patmatmotifs
   loop will be slow as it processes one sequence at a time.
   This could be parallelised with concurrent.futures if needed.

5. The pipeline has no resume capability. If it fails midway,
   it must be restarted from the beginning. Intermediate files
   in the output directory are not reused.

================================================================
