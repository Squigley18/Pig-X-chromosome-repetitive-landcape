#!/bin/bash
# File: PlotWrapper.sh
#
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #Call functions from shell-script-functions.sh script 

parameter_check $1
cd $1 

create_folder overlap.hits.to.chromosome  
create_folder segments
create_folder ideograms
create_folder tables 

/home/sq16564/Analysis-of-pig-X/./overlap.hits.entire.chromosome.R identity_99 identity_95 repeatmasked_99 repeatmasked_95  #Run overlap.hits.entire.chromosome.R script with output from LASTZ alignment tables 

mv *.png overlap.hits.to.chromosome

/home/sq16564/Analysis-of-pig-X/./entire.chromosome.R identity_99 identity_95 repeatmasked_99 repeatmasked_95 #Run entire.chromosome.R script with output from LASTZ alignment tables

mv *ideogram.png ideograms
mv *segments.svg segments

for t in *.tab; do perl /home/sq16564/Analysis-of-pig-X/./perl.extract.tsv.files.pl $t; done #for each table produced by the entire.chromosome.R with the hit start and stop positions run through the perl.extract.tsv.files.pl script to extract the FASTA sequence within those regions

mv *.tab tables

#Run the DNA sequences from the perl script through BLAST to provide a table of subject sequences which then can be used to calculate the frequency of subject sequences, and filtering out the low complexity DNA, and determine the location of the subject sequence in the X chromosome and extract the sequence for manual BLASTING
/home/sq16564/Analysis-of-pig-X/./full.blast.filtering.sh Rm99 100 
/home/sq16564/Analysis-of-pig-X/./full.blast.filtering.sh ID99 1000
/home/sq16564/Analysis-of-pig-X/./full.blast.filtering.sh Rm95 2000
/home/sq16564/Analysis-of-pig-X/./full.blast.filtering.sh ID95 10000
