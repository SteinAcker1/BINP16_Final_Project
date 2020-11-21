Title: restriction_app.py
Author: Stein Acker
Date: 28-10-2020

***IMPORTANT***
For this app to work properly, the files restriction_app.py,
restriction_sites.tsv, and project_module.py MUST be in the same directory. To
use the optional helper interface rather than the argparse interface,
restriction_helper.sh must also be in the directory.

Description:
    This program will, given a FASTA or FASTQ file, generate a list of
    potential restriction sites in the file as well as their corresponding
    positions. Optional user inputs include a specific set of enzymes to
    search and a specific set of substrings to filter for in headers.

List of user-defined functions:
    get_complement(seq): Produces the complementary sequence of DNA in the
    5'-3' direction

    create_DNA_regex(seq): Translates IUPAC nucleotide codes into
    computer-readable regular expressions
    get_regex_complement(seq): Gets the complementary DNA sequence and
    converts IUPAC codes to regex in one command

    TSV_to_list(file): Converts a TSV file to a Python list, then converts
    IUPAC to regex and eliminates all non-nucleotide or regex characters from
    a DNA sequence

    add_complements_to_list(ls): Adds the complements of all non-palindromic
    restriction sites to the list.

    fastq_to_fasta(fastq): Converts FASTQ files to FASTA format

    single_line_fasta(fasta): Converts multiline FASTA files to single-line
    FASTA files

List of non-standard modules:
    project_module: Contains all of the custom functions used in this app

Procedure:
    Step 1: Load inputs and format them correctly
    Step 2: If the input is in FASTQ format, convert to FASTA
    Step 3: Filtering only for inputted headers and enzymes, iterate through
    the lines of the FASTA file and generate a list of restriction site
    matches using re.finditer()
    Step 4: Print a list of found restriction sites, their positions, their
    enzymes, and their illustrated sites in an easy-to-read fashion

Usage:

Helper (recommended if not comfortable with argparse):
  $bash restriction_helper.sh
  Follow the prompts as they appear in the terminal

argparse:
    $python3 restriction_app.py
    -i INFILE

  -k KEYWORD            Allows the user to enter key phrases to search for in
                        the FASTA/FASTQ headers. Leave blank to include all
                        headers. Key phrases must be comma-separated with no
                        spaces. To search for a specific header, enter in the
                        format exact,headerName. Case insensitive.

  -e ENZYMES            List of restriction enzymes to filter for.
                        Case insensitive. Separate with commas, but do not
                        include spaces. Partial name also works if you want to
                        return multiple enzymes (e.g. eco for EcoRI and EcoRV).
                        Enter "all" or leave blank to search restriction sites
                        for all enzymes.

  -r RESTRICTION SITES TSV
                        Allows the user to enter the filepath to their own TSV
                        file containing restriction site sequences and
                        restriction enzymes. If left blank, the default TSV will
                        be used. Look at the default TSV to ensure yours is
                        formatted correctly.
  -o OUTFILE



  Source for default restriction sites TSV:
  https://international.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities
