#!/bin/bash

#$ -cwd
#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #Call the shell script functions pre-prepared in the shell-script-functions.sh script

create_folder tsv_files
create_folder blast.tables
create_folder refined.blasts

for f in *.out.tsv #all files from perl extraction of FASTA DNA sequence in the palindrome loocation on the X chromosome
do
 log "blasting $f";
 blastn -db /storage/projects/SkinnerLab/nt/nt -query $f  -num_threads 8 -out $f.blast -outfmt "6 qseqid sseqid  pident sstart send qstart qend stitle sacc length evalue"
#-db where the blast database is stored, -query the DNA sequences to be blasted, -num_threads computing power needed,
# -out the output file written to, -outfmt is the output format where 6 produces tabular format, qseqid query sequence id,
# sseqid subject sequence id, pident percentage identity,sstart subjects start, send subject end, qstart query start,
 #qend query end, stitle subject title, sacc subject accession, length alignment length, evalue expect value
done

 for f in *.blast; do awk '!/clone/' $f > $f.out ; done #remove all unwanted alignments to clone sequences

 for f in *.blast.out; do  awk '!/breed/' $f > $f.final; done #remove all unwanted alignments to breed sequences

 rm *.out #remove first output file with only clone sequences removed

 for f in *.out.tsv.*; do mv "$f" "${f/.out.tsv./.}"; done #rename output files to remove .out.tsv from name

mv *.final refined.blasts

mv *.out.tsv tsv_files
mv *.blast blast.tables

create_folder unmasked.blast
create_folder masked.blast

for f in refined.blasts/all.unmasked.**.final; do cut -f9 $f >> unmasked.accessions; done #create a file of all the accession numbers for the unmasked palindrome subject sequences found in the blast search
for f in refined.blasts/all.masked.**.final; do cut -f9 $f >> masked.accessions; done #create a file of all the accession numbers for the masked palindrome subject sequences found in the blast search

/home/sq16564/Analysis-of-pig-X/./frequency.tables.R unmasked.accessions #Calculate the frequency of occurances for each unique accession number for the unmasked palindrome blast search subject sequences
/home/sq16564/Analysis-of-pig-X/./frequency.tables.R masked.accessions #Calculate the frequency of occurances for each unique accession number for the masked palindrome blast search subject sequences

/home/sq16564/Analysis-of-pig-X/./frequency.bar.graph.R unmasked.accessions.freq #Create a bar graph showing the frequency of occurances for each unique accession number in the unmasked palindrome blast search subject sequences
create_folder unmasked.blast/bar.graphs
mv *.png unmasked.blast/bar.graphs
/home/sq16564/Analysis-of-pig-X/./frequency.bar.graph.R masked.accessions.freq #Create a bar graph showing the frequency of occurances for each unique accession number in the masked palindrome blast search subject sequences
create_folder masked.blast/bar.graphs
mv *.png masked.blast/bar.graphs

awk '(NR==1) || ($2 > 2 ) ' unmasked.accessions.freq > highest.freq.unmasked #Filter the accession numbers from the unmasked palindrome blast search subject sequences to remove any hits with below 2 occurances
awk '(NR==1) || ($2 > 2 ) ' masked.accessions.freq > highest.freq.masked #Filter the accession numbers from the masked palindrome blast search subject sequences to remove any hits with below 2 occurances

cut -f1 highest.freq.unmasked > highest.accessions.unmasked #Isolate the highest frequency accession numbers from the unmasked palindrome blast search subject sequences
cut -f1 highest.freq.masked > highest.accessions.masked #Isolate the highest frequency accession numbers from the masked palindrome blast search subject sequences

rm highest.freq.unmasked
rm highest.freq.masked

cat refined.blasts/unmasked**.final >> unmasked.combined.blast #Make file with combined blast results for unmasked palindromes for use later in analysis
cat refined.blasts/masked**.final >> masked.combined.blast #Make file with combined blast results for masked palindromes for use later in analysis

cut -f8,9 masked.combined.blast > masked.combined.blast.titles #extract the accession numbers and subject titles from the combined blast results of the unmasked palindromes
cut -f8,9 unmasked.combined.blast > unmasked.combined.blast.titles #extract the accession numbers and subject titles from the combined blast results of the masked palindromes

rm masked.combined.blast
rm unmasked.combined.blast

awk  '{ if (!a[$3]++ ) print ;}' unmasked.combined.blast.titles > unmasked.unique.blast.titles  #remove duplicates of unique accession sequences from unmasked palindrome blast search subject sequences
awk  '{ if (!a[$3]++ ) print ;}' masked.combined.blast.titles > masked.unique.blast.titles #remove duplicates of unique accession sequences from masked palindrome blast search subject sequences

rm unmasked.combined.blast.titles 
rm masked.combined.blast.titles

grep -f highest.accessions.unmasked unmasked.unique.blast.titles > unmasked.most.frequent.blast.titles #Using the highest frequency accession numbers to extract the titles for the most frequent blast results 
grep -f highest.accessions.masked masked.unique.blast.titles > masked.most.frequent.blast.titles #Using the highest frequency accession numbers to extract the titles for the most frequent blast results 

create_folder unmasked.blast/subject.tables
create_folder unmasked.blast/subject.chromosome.positions
create_folder unmasked.blast/subject.tsv.files
create_folder unmasked.blast/subjects

while read -r line #For every accession number in the highest frequency accessions do the following:
do
for f in refined.blasts/all.unmasked.**.final; do grep $line $f >> $line; done #extract the lines with the corresponding accession number from the blast tables of results
cut -f1,4,5,6,7 $line > $line.position #extract the palindrome hit data, query start and stop position, and subject start and stop position
sed 's/:/\t/g' $line.position > tab #replace : divisions with tabs to allow for extraction of required data
cut -f5,6,8,9,10,11 tab > $line.tab #Take the palindrome start and stop position, query start and stop position, and subject start and stop position
/home/sq16564/Analysis-of-pig-X/./blast.subject.tables.R $line.tab #Calculate where in the X chromosome the homology to the subject sequence is found
rm $line.position
rm tab
split -l 100 $line.tab.chromosome.hits $line.chromosome. #split the tables from blast.subject.tables.R script to ensure they are not too large for manual blasting when extracting FASTA DNA sequences
for l in $line.chromosome.*; do perl /home/sq16564/Analysis-of-pig-X/./perl.extract.tsv.files.pl $l;done #Using the perl script extract the FASTA DNA files from the pig X chromosome from within the regions of homology to the blast subject sequences
mv $line unmasked.blast/subjects
mv $line.tab unmasked.blast/subject.tables
mv $line.tab.blast.hits unmasked.blast/subject.tables
mv $line.tab.chromosome.hits unmasked.blast/subject.chromosome.positions
mv *.out.tsv unmasked.blast/subject.tsv.files
rm $line.chromosome.*
done < highest.accessions.unmasked

create_folder masked.blast/subject.tables
create_folder masked.blast/subject.chromosome.positions
create_folder masked.blast/subject.tsv.files
create_folder masked.blast/subjects

while read -r line #For every accession number in the highest frequency accessions do the following:
do
for f in refined.blasts/all.masked.**.final; do grep $line $f >> $line; done #extract the lines with the corresponding accession number from the blast tables of results
cut -f1,4,5,6,7 $line > $line.position  #extract the palindrome hit data, query start and stop position, and subject start and stop position
sed 's/:/\t/g' $line.position > tab #replace : divisions with tabs to allow for extraction of required data
cut -f5,6,8,9,10,11 tab > $line.tab #Take the palindrome start and stop position, query start and stop position, and subject start and stop position
/home/sq16564/Analysis-of-pig-X/./blast.subject.tables.R $line.tab #Calculate where in the X chromosome the homology to the subject sequence is found
rm $line.position
rm tab
split -l 100 $line.tab.chromosome.hits $line.chromosome. #split the tables from blast.subject.tables.R script to ensure they are not too large for manual blasting when extracting FASTA DNA sequences
for l in $line.chromosome.*; do perl /home/sq16564/Analysis-of-pig-X/./perl.extract.tsv.files.pl $l;done #Using the perl script extract the FASTA DNA files from the pig X chromosome from within the regions of homology to the blast subject sequences
mv $line masked.blast/subjects
mv $line.tab masked.blast/subject.tables
mv $line.tab.blast.hits masked.blast/subject.tables
mv $line.tab.chromosome.hits masked.blast/subject.chromosome.positions
mv *.out.tsv masked.blast/subject.tsv.files
rm $line.chromosome.*
done < highest.accessions.masked
