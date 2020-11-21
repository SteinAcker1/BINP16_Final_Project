#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 10:53:39 2020

@author: steinacker
"""
import re

def get_complement(seq):
    #Produces the complementary sequence of DNA in the 5'-3' direction
    table = seq.maketrans('ACTGactg[]','TGACtgac][')
    comp = seq.translate(table)[::-1]
    return comp

def create_DNA_regex(seq):
    #Translates IUPAC nucleotide codes into computer-readable regular expressions
    seq = seq.upper()
    symbols = {ord('Y'):'[CT]',
               ord('R'):'[AG]',
               ord('S'):'[CG]',
               ord('W'):'[AT]',
               ord('K'):'[GT]',
               ord('M'):'[CA]',
               ord('B'):'[CGT]',
               ord('D'):'[AGT]',
               ord('H'):'[CAT]',
               ord('V'):'[CGA]',
               ord('N'):'[ACGT]'}
    return seq.translate(symbols)

def get_regex_complement(seq):
    #Gets the complementary DNA sequence and converts IUPAC codes to regex in one command
    return get_complement(create_DNA_regex(seq))

def tidy_up(seq):
    #Converts IUPAC to regex and eliminates all non-nucleotide or regex characters from a DNA sequence
        seq = create_DNA_regex(seq)
        seq = re.sub('[^ACTGactg\[\]]','',seq)
        return seq

def TSV_to_list(file):
    #Converts a TSV file to a Python list, then converts IUPAC to regex and eliminates all non-nucleotide or regex characters from a DNA sequence
    ls = []
    for line in file:
        line = line.strip('\n')
        line = line.split('\t')
        line[0] = tidy_up(line[0])
        if len(line) > 1:
            ls.append(line)
    return ls

def convert_list_to_dict(ls):
    #Converts a list to a dictionary. Not used in project.
    dictionary = dict()
    for line in ls:
        dictionary[line[0]] = line[1]
    return dictionary

def add_complements_to_list(ls):
    #Adds the complements of all non-palindromic restriction sites to the list.
    seqs = []
    for line in ls:
        seqs.append(line[0])
    for line in ls:
        comp = get_complement(line[0])
        if comp not in seqs:
            ls.append([comp,line[1]])
    return ls

def fastq_to_fasta(fastq):
    #Converts FASTQ files to FASTA format
    fasta = []
    for line in range(len(fastq)):
        if fastq[line][0] == '@' and line != len(fastq)-1:
            if not re.search('[^ACTGactg]',fastq[line+1].strip('\n')):
                fasta.append('>' + fastq[line][1:])
        elif not re.search('[^ACTGactg]',fastq[line].strip('\n')):
            fasta.append(fastq[line])
    return fasta

def single_line_fasta(fasta):
    #Converts multiline FASTA files to single-line FASTA files
    newFasta = []
    seq = ''
    for line in fasta:
        if line[0] == '>':
            if len(newFasta) > 0:
                newFasta.append(seq)
            newFasta.append(line)
            seq = ''
        elif len(seq) > 0:
            line = line.strip('\n')
            seq += line
        else:
            seq = line.strip('\n')
    newFasta.append(seq)
    return newFasta
