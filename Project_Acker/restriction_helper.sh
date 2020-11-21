#!/bin/sh

#Part 1: Introducing prompts and variables
echo Welcome to Restriction App! Please enter the filepath to the FASTA/FASTQ
echo file you would like to analyze:
read infile
echo \-----

echo Enter the restriction enzymes you would like to search for, using only
echo commas as seperators. Case insensitive and partial names are fine \(e.g.
echo eco,m,k for enzymes containing \"eco\", \"m\", or \"k\"\). Leave blank to
echo search for all enzymes:
read enzymes
echo \-----

echo Enter the key phrases you would like to search for in the headers. Same
echo rules as above. To include only exact matches, begin with \'exact\'
echo and a comma \(e.g. exact,scaffold_5 will only search the header named
echo scaffold_5, while scaffold_5 will also search scaffold_50, etc.\). Leave
echo blank to search for all headers:
read phrases
echo \-----

echo Enter a custom TSV file for restriction enzymes or leave blank for the
echo default TSV file:
read TSV
echo \-----

echo Finally, enter a filepath for your output file or leave blank to print the
echo results to the terminal:
read output
echo \-----

#Part 2: Converting variables to parsable arguments
infile=$(echo \-i $infile)

if [[ $enzymes != '' ]]
then
  enzymes=$(echo \-z $enzymes)
else
  enzymes=$(echo '')
fi

if [[ $phrases != '' ]]
then
  phrases=$(echo \-k $phrases)
else
  phrases=$(echo '')
fi

if [[ $TSV != '' ]]
then
  TSV=$(echo \-r $TSV)
else
  TSV=$(echo '')
fi

if [[ $output != '' ]]
then
  output=$(echo \-o $output)
else
  output=$(echo '')
fi

#Part 3: Executing the program
echo Loading...
eval python3 restriction_app.py $infile $enzymes $phrases $TSV $output
echo Done!
