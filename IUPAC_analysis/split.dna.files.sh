#!/bin/bash
# File: PlotWrapper.sh
#
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #Call the the shell script functions pre-defined in the script

#The pig DNA sequence file is around 130MB in size therefore splittling the file into two 65MB files will be required to reduce the computational demand when performing the IUPACpal software 
split --bytes=65MB ~/dna_files/Sus_scrofa.Sscrofa11.1.dna.chromosome.X.fa pig. 
split --bytes=65MB ~/dna_files/Sus_scrofa.Sscrofa11.1.dna_rm.chromosome.X.fa pig.rm.

#Remove the FASTA header from the files 
tail -n +2 pig.aa > pig
rm pig.aa
tail -n +2 pig.rm.aa > pig.rm
rm pig.rm.aa

#Add a new header in FASTA style to each split DNA file which will be required for the IUPACpal software analysis
echo ">X dna:chromosome:Sscrofa11.1:first.half" > Sus_scrofa11.1.dna.X.first.half.fa
echo ">X dna:chromosome:Sscrofa11.1:second.half" > Sus_scrofa11.1.dna.X.second.half.fa
echo ">X dna:chromosome:Sscrofa11.1:rm:first.half" > Sus_scrofa11.1.dna.X.rm.first.half.fa
echo ">X dna:chromosome:Sscrofa11.1:rm:second.half" > Sus_scrofa11.1.dna.X.rm.second.half.fa

cat pig >> Sus_scrofa11.1.dna.X.first.half.fa
cat pig.rm >> Sus_scrofa11.1.dna.X.rm.first.half.fa
cat pig.ab >> Sus_scrofa11.1.dna.X.second.half.fa
cat pig.rm.ab >> Sus_scrofa11.1.dna.X.rm.second.half.fa

rm pig
rm pig.rm
rm pig.ab
rm pig.rm.ab
