#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: restriction_app.py
Author: Stein Acker
Date: 28-10-2020

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
    $python3 restriction_app.py
    -i INFILE

  -k KEYWORD            Allows the user to enter key phrases to search for in
                        the FASTA/FASTQ headers. Leave blank to include all
                        headers. Key phrases must be comma-separated with no
                        spaces. To search for a specific header, enter in the
                        format exact,headerName. Case insensitive.

  -z ENZYMES            List of restriction enzymes to filter for.
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

"""
#Setup
import argparse
import sys
import re
import copy as c
import project_module as p

usage = 'This program will, given a FASTA or FASTQ file, generate a list of potential restriction sites in the file as well as their corresponding positions.'
parser = argparse.ArgumentParser(description = usage)
parser.add_argument('-i',
                    dest = 'infile',
                    metavar = 'INFILE',
                    type = argparse.FileType('r'),
                    required = True
                    )
parser.add_argument('-k',
                    dest = 'keyword',
                    metavar = 'KEYWORD',
                    type = str,
                    help = 'Allows the user to enter key phrases to search for in the FASTA/FASTQ headers. Leave blank to include all headers. Key phrases must be comma-separated with no spaces. To search for a specific header, enter in the format exact,headerName. Case insensitive.',
                    default = '>'
                    )
parser.add_argument('-z',
                    dest = 'enzymes',
                    metavar = 'ENZYMES',
                    type = str,
                    help = 'List of restriction enzymes to filter for. Case insensitive. Separate with commas, but do not include spaces. Partial name also works if you want to return multiple enzymes (e.g. eco for EcoRI and EcoRV). Enter "all" or leave blank to search restriction sites for all enzymes.',
                    default = 'all'
                    )
parser.add_argument('-r',
                    dest = 'rsites',
                    metavar = 'RESTRICTION SITES TSV',
                    type = str,
                    help = 'Allows the user to enter the filepath to their own TSV file containing restriction site sequences and restriction enzymes. If left blank, the default TSV will be used. Look at the default TSV to ensure yours is formatted correctly.',
                    default = 'restriction_sites.tsv'
                    )
parser.add_argument('-o',
                    dest = 'outfile',
                    metavar = 'OUTFILE',
                    type = argparse.FileType('w'),
                    default = sys.stdout
                    )

args = parser.parse_args()

#Opening up the TSV
restrict_data = args.rsites

with open(args.rsites, 'r') as tsv:
    all_restricts = p.TSV_to_list(tsv)
    all_restricts = p.add_complements_to_list(all_restricts)

#Reading the enzyme input to decide which enzymes to look for
enzyme_input = args.enzymes.upper().split(',')
#If the user asks for all enzymes, the "restricts" variable is simply a deepcopy of all_restricts
if enzyme_input[0] == 'ALL':
    restricts = c.deepcopy(all_restricts)
else:
    restricts = []
    for line in all_restricts:
        for enzyme in enzyme_input:
            if enzyme in line[1].upper():
                restricts.append(line)

#Deleting all_restricts to clear up memory
del all_restricts

#Open the input and output
with args.infile as i, args.outfile as o:
    #Convert the input to a list
    file = i.readlines()
    #If the input is in FASTQ format, then convert it to FASTA format
    if file[0][0] == '@':
        file = p.fastq_to_fasta(file)
    #Convert the FASTA file to single-line format
    file = p.single_line_fasta(file)
    #Set some marker variables
    flag = False
    sites = False
    #Turn the keyword input into a list to iterate through
    headers = args.keyword.upper().split(',')
    #For each line in the FASTA...
    for line in file:
        #...and for each header mentioned (default being >, or all headers)...
        for header in headers:
            #if the line is a header and the keyword is in the line OR the header is equivalent to the keyword...
            if line[0] == '>':
                #save that line as "ID" and set flag to equal True
                if headers[0] == 'EXACT' and line.upper().strip('>').strip('\n') == header:
                    ID = line[1:].strip('\n')
                    flag = True
                elif headers[0] != 'EXACT' and header in line.upper():
                    ID = line[1:].strip('\n')
                    flag = True
            #If flag is true...
            elif flag:
                #start a repository for sites already used
                used_sites = []
                #For each enzyme being tested...
                for enzyme in restricts:
                    #for each restriction site match found in the line...
                    for match in re.finditer(enzyme[0],line):
                        #define the start and end positions of the restriction sites, as well as the site and its complement...
                        start = match.start()
                        end = match.end()
                        seq = line[start:end]
                        comp = p.get_complement(seq)[::-1]
                        #and create an ID for this site.
                        string = str(start) + enzyme[1] + str(end)
                        #If the ID is not in the used sites repository...
                        if string not in used_sites:
                            #add this entry to the output and add the used ID to the used sites repository.
                            o.write('---------------------\n\n' +
                                'Sequence: {}\n'.format(ID) +
                                'Enzyme: {}\n'.format(enzyme[1]) +
                                'Position: {}-{}\n\n'.format(start + 1,end + 1) +
                                '5\'...{}...3\'\n'.format(seq) +
                                '3\'...{}...5\'\n\n'.format(comp) +
                                '---------------------\n')
                            used_sites.append(string)
                        #Set sites to True to signal that restriction sites were found
                        sites = True
                        #Reset the flag variable to False
                        flag = False
    #If the sites variable remains false, output "No restriction sites were found"
    if sites == False:
        o.write('No restriction sites found\n')
