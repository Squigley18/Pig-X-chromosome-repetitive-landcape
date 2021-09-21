#!/bin/bash
# File: PlotWrapper.sh
#
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #Call shell script functions as defined in the script

./split.dna.files.sh #Split the DNA files in two to reduce the computational demand required to run the software

./iupacpal.sh Sus_scrofa11.1.dna.X.rm.first.half.fa "X dna:chromosome:Sscrofa11.1:rm:first.half" masked.1 #Run IUPACpal on the first half of the masked DNA sequence
./iupacpal.sh Sus_scrofa11.1.dna.X.rm.second.half.fa "X dna:chromosome:Sscrofa11.1:rm:second.half" masked.2 #Run IUPACpal on the second half of the masked DNA sequence
./iupacpal.sh Sus_scrofa11.1.dna.X.first.half.fa "X dna:chromosome:Sscrofa11.1:first.half" unmasked.1 #Run IUPACpal on the first half of the unmasked DNA sequence
./iupacpal.sh Sus_scrofa11.1.dna.X.second.half.fa "X dna:chromosome:Sscrofa11.1:second.half" unmasked.2 #Run IUPACpal on the second half of the unmasked DNA sequence

create_folder dna_files
mv *.fa dna_files

./iupac.table.sh unmasked.1.iupac unmasked.2.iupac #Convert the text format of the iupac output into a table which can be maniputlated for further analysis
./iupac.table.masked.sh masked.1.iupac masked.2.iupac #Convert the text format of the iupac output into a table which can be maniputlated for further analysis

create_folder iupac.output
mv *.iupac iupac.output

create_folder chromosome.overlaps

./overlap.hits.entire.chr.iupac.R unmasked.1.iupac.tab masked.2.iupac.tab  #Create a graphic representation of where the palindromes are found on the chromosome

mv *.png chromosome.overlaps

create_folder segments
create_folder ideograms
create_folder tables

./entire.chromosome.iupac.R  unmasked.1.iupac.tab masked.2.iupac.tab #Create a graphic representation of the genes the palindromes overlap, where the complimentary palindromes are found and export a table of palindrome start and sttop postitions for extracting the DNA FASTA format

mv *segments.svg segments
mv *ideogram.png ideograms

for t in *.pal.tab; do perl /home/sq16564/Analysis-of-pig-X/./perl.extract.tsv.files.pl $t; done #Using the perl script and the table of palindrome positions the X chromosome FASTA sequence was extracted to be run through BLAST

mv *.pal.tab tables 

./iupac.blast.filtering.sh #Running the FASTA DNA sequence through BLAST returns a table of subject sequences showing homology to the DNA sequence then the exact region of the homology is isolated for further manual BLASTING 
