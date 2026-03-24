================================================================
  Protein Family Conservation & Motif Pipeline
  USER MANUAL — Section 1: Guide for the Non-Programmer
================================================================

WHAT DOES THIS PROGRAMME DO?
------------------------------
This programme is an automated tool for biologists. You give it
the name of a protein family and a group of organisms, and it
will do the following for you automatically:

  1. Download matching protein sequences from the NCBI database
     (the same database used by BLAST on the NCBI website).

  2. Tell you how many sequences were found and which species
     they come from, so YOU can decide whether it is worth
     continuing with that dataset.

  3. Align all the sequences and produce a CONSERVATION PLOT
     showing which parts of the protein are most similar across
     different species.

  4. Scan every sequence for known protein MOTIFS (functional
     domains) from the PROSITE database, and report which motifs
     were found and how often.

  5. Optionally calculate physicochemical properties of your
     protein sequences (molecular weight, isoelectric point).

  6. Save all results neatly into a folder you name yourself.

----------------------------------------------------------------
BEFORE YOU START
----------------------------------------------------------------
You must run this programme on the MSc bioinformatics server,
where all the required tools are already installed.

You will need:
  - A terminal window logged in to the MSc server
  - The file protein_pipeline.py saved in your home directory
  - An internet connection (the server uses this to contact NCBI)

You do NOT need to install anything yourself.

----------------------------------------------------------------
HOW TO RUN THE PROGRAMME
----------------------------------------------------------------
1. Open a terminal and log in to the MSc server.

2. Check that protein_pipeline.py is in your current directory:

       ls

   You should see protein_pipeline.py listed.

3. Start the programme by typing:

       python3 protein_pipeline.py

   Then press Enter.

4. The programme will greet you with a short description and
   then ask you three questions. Answer each one and press
   Enter to continue.

----------------------------------------------------------------
QUESTION 1: Protein family name
----------------------------------------------------------------
Type the name of the protein family you want to study.

Good examples:
  glucose-6-phosphatase
  ABC transporter
  adenylyl cyclase
  protein kinase C

Tips:
  - Use the common scientific name, not a gene symbol
  - If you get no results, try a shorter or simpler name
    (e.g. "kinase" instead of "serine/threonine protein kinase")
  - Spelling must be correct - the programme searches NCBI
    exactly as you type it

*** TEST EXAMPLE: glucose-6-phosphatase ***

----------------------------------------------------------------
QUESTION 2: Taxonomic group
----------------------------------------------------------------
Type the scientific name of the group of organisms you want to
search within.

Good examples:
  Aves           (all birds)
  Mammalia       (all mammals)
  Rodentia       (rodents)
  Vertebrata     (all vertebrates)
  Insecta        (insects)
  Actinopterygii (ray-finned fish)

Tips:
  - Use the Latin/scientific name, not the common name
    ("Aves" works; "birds" may not)
  - If you get too many sequences (more than 1,000), try a
    smaller group (e.g. "Galliformes" instead of "Aves")
  - If you get too few sequences (fewer than 5 species), try
    a larger group (e.g. "Vertebrata" instead of "Aves")

*** TEST EXAMPLE: Aves ***

----------------------------------------------------------------
QUESTION 3: Output directory name
----------------------------------------------------------------
Type a name for the folder where all your results will be saved.
If you leave this blank and just press Enter, the folder will
be called "pipeline_output".

Good examples:
  g6pase_birds
  kinase_mammals_results
  my_project_output

Tips:
  - Use underscores instead of spaces in the folder name
  - If a folder with that name already exists, the programme
    will warn you and ask whether to continue

*** TEST EXAMPLE: g6pase_birds ***

----------------------------------------------------------------
WHAT HAPPENS WHILE THE PROGRAMME RUNS
----------------------------------------------------------------
After you answer the three questions, the programme runs
automatically. You will see status messages like:

  [Running] esearch -db protein ...
  [OK] Sequences saved to: g6pase_birds/sequences_raw.fasta

These messages tell you what the programme is doing at each
stage. You do not need to do anything - just wait.

At ONE point during the run, the programme will PAUSE and ask
you a question:

  PAUSE 1 - After downloading sequences:
  ---------------------------------------
  The programme will show you a table of how many sequences
  come from each species. For example:

      Total sequences : 47
      Unique species  : 19

      Species                                       Count
      -----------------------------------------------------
      Gallus gallus                                    12
      Taeniopygia guttata                               6
      Columba livia                                     4
      Meleagris gallopavo                               3
      ...

  It will then ask:

      Continue with this dataset? [y/N]:

  Type  y  and press Enter to continue.
  Type  n  and press Enter to stop and try a different search.

  If only 1 or 2 species are present, the programme will warn
  you - conservation analysis works best with many species.

  PAUSE 2 - After the main analysis:
  ------------------------------------
  The programme will ask:

      Run optional EMBOSS pepstats analysis? [y/N]:

  Type  y  if you want molecular weight and isoelectric point
  data for your sequences. Type  n  to skip this step.
  Either choice is fine.

----------------------------------------------------------------
YOUR RESULTS - WHAT EACH FILE CONTAINS
----------------------------------------------------------------
All results are saved in the folder you named. Here is what
each output file contains and what it means:

  pipeline_summary.txt
  ---------------------
  A brief record of the whole run: the query used, how many
  sequences were found, and a list of all output files.
  READ THIS FIRST if you want a quick overview.

  species_diversity.txt
  ----------------------
  A complete list of all species found in your dataset, with
  the number of sequences from each species.
  Useful for checking whether your dataset covers many species
  or is dominated by one or two.

  sequences_raw.fasta
  --------------------
  All the protein sequences downloaded from NCBI, in FASTA
  format. Each entry starts with a ">" line containing the
  sequence name and species. You can open this in a text editor
  or upload it to NCBI BLAST if you wish.

  sequences_for_alignment.fasta
  ------------------------------
  The sequences that were actually used for the alignment step.
  If more than 100 sequences were downloaded, the programme
  selected a representative set of 100 (choosing the longest,
  most complete sequences). Otherwise, this is the same as
  sequences_raw.fasta.

  sequences_aligned.fasta
  ------------------------
  The aligned sequences. Dashes (-) have been inserted to line
  up matching regions across all sequences. This is what the
  conservation analysis is based on.
  You do not normally need to look at this file directly.

  conservation_plot.png
  ----------------------
  *** THIS IS ONE OF THE MOST IMPORTANT OUTPUT FILES ***

  A graph (image file) showing sequence conservation along the
  full length of the protein. The x-axis is the position along
  the protein sequence. The y-axis is the conservation score.

  HOW TO INTERPRET IT:
    - HIGH peaks = positions where the amino acid is the same
      (or very similar) across most species. These positions
      are likely to be functionally critical - evolution has
      kept them the same because changing them would break the
      protein.
    - LOW troughs = positions where the amino acid varies
      between species. These regions may be less critical, or
      may be involved in species-specific adaptations.
    - A protein with mostly HIGH conservation across Aves is
      likely performing the same essential function in all
      birds.

  Open this file with any image viewer (e.g. type:
      display conservation_plot.png
  or copy it to your local computer).

  conservation_scores.txt
  ------------------------
  The numbers behind the conservation plot, saved as a table.
  Column 1 = position number in the alignment.
  Column 2 = conservation score (0.00 = fully variable,
                                  1.00 = identical in all seqs).
  You can import this into Excel or R for further plotting
  if you wish.

  prosite_summary.txt
  --------------------
  *** THIS IS THE OTHER MOST IMPORTANT OUTPUT FILE ***

  A table of PROSITE motifs (known functional domains) found
  in your sequences. Each row shows:
    - The PROSITE motif ID and name (e.g. PS00174 G6PASE)
    - How many sequences it was found in

  HOW TO INTERPRET IT:
    - If a motif appears in MOST or ALL of your sequences,
      it is likely a defining structural/functional feature
      of this protein family.
    - If a motif appears in only a few sequences, it may be
      a secondary feature or present only in certain species.
    - If NO motifs are found, try a broader search, or the
      protein family may simply not be well-represented in
      the PROSITE database yet.

  prosite_results/  (folder)
  ---------------------------
  A folder containing the individual PROSITE scan result for
  every sequence. You normally do not need to look at these
  directly - the summary file above contains the key findings.

  pepstats_results.txt  (only if you chose to run pepstats)
  -----------------------------------------------------------
  Physicochemical properties for each sequence:
    - Molecular weight (in Daltons, Da)
    - Isoelectric point (pI) - the pH at which the protein
      carries no net charge
    - Amino acid composition (percentage of each amino acid)

  The programme also prints a summary of mean, minimum, and
  maximum MW and pI values to the screen.

----------------------------------------------------------------
EXAMPLE: glucose-6-phosphatase in Aves (birds)
----------------------------------------------------------------
Input given to the programme:
  Protein family name : glucose-6-phosphatase
  Taxonomic group     : Aves
  Output directory    : g6pase_birds

What to expect:
  - Approximately 20-60 sequences retrieved from NCBI
    (the exact number depends on current NCBI content)
  - Multiple bird species represented, e.g.:
      Gallus gallus (chicken)
      Taeniopygia guttata (zebra finch)
      Columba livia (pigeon)
      Meleagris gallopavo (turkey)
      Anas platyrhynchos (duck)
  - Conservation plot shows HIGH conservation in the catalytic
    domain region of the protein, reflecting its essential
    metabolic role (glucose-6-phosphatase is critical for
    blood glucose regulation in all vertebrates)
  - PROSITE scan finds the glucose-6-phosphatase active site
    motif (PS00174) in most or all sequences
  - pepstats shows molecular weight in the ~36-40 kDa range
    (consistent with the known size of this enzyme)

----------------------------------------------------------------
COMMON PROBLEMS AND SOLUTIONS
----------------------------------------------------------------

PROBLEM: "No sequences found for this query."
  CAUSE:   The protein name or taxon was not recognised by NCBI,
           or there are genuinely no sequences for that combination.
  SOLUTION:
    - Double-check the spelling of both the protein name and taxon
    - Try a shorter protein name (e.g. "phosphatase" instead of
      "glucose-6-phosphatase catalytic subunit alpha")
    - Try a broader taxon (e.g. "Vertebrata" instead of "Aves")
    - Search for the protein manually on the NCBI website
      (https://www.ncbi.nlm.nih.gov/protein/) to check what
      name NCBI uses

PROBLEM: "Program 'esearch' not found on PATH."
  CAUSE:   A required tool is not installed or not available.
  SOLUTION: Contact your system administrator. This programme
            must be run on the MSc server where tools are installed.

PROBLEM: Warning that only 1 or 2 species are represented.
  CAUSE:   The chosen taxon is very narrow, or the protein is
           only sequenced in a few species.
  SOLUTION: Use a broader taxonomic group for more species.
            Conservation analysis is most informative with
            sequences from 10 or more species.

PROBLEM: More than 1,000 sequences found.
  CAUSE:   Very common protein family or very large taxon.
  SOLUTION: The programme will ask if you want to continue
            with the first 1,000 sequences. You can type y,
            or type n and re-run with a more specific protein
            name or a narrower taxonomic group.

PROBLEM: The output folder already exists.
  CAUSE:   You have run the programme before with the same
           output directory name.
  SOLUTION: Either type y to overwrite the previous results,
            or choose a new directory name to keep both sets
            of results.

PROBLEM: No PROSITE motifs found.
  CAUSE:   The protein family may not be in the PROSITE
           database, or the sequences may be too divergent
           for pattern matching.
  SOLUTION: This is a valid biological result - report it as
            "no known PROSITE motifs identified". You can
            search the PROSITE database manually at
            https://prosite.expasy.org/ for further information.

----------------------------------------------------------------
TIPS FOR GETTING THE BEST RESULTS
----------------------------------------------------------------
  - Use scientific names for the taxonomic group (Aves, not birds)
  - Aim for a dataset with 10-100 sequences from 5+ species
  - If you get more than 500 sequences, the analysis may take
    a long time - consider narrowing the search
  - The conservation plot is most informative when sequences
    come from DIVERSE species (e.g. different bird orders),
    not just multiple sequences from one species
  - Save your results folder somewhere safe before re-running,
    as re-running to the same directory will overwrite files

================================================================
