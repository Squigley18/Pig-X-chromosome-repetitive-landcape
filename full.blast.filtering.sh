#!/bin/bash

#$ -cwd
#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #call the script containing the shell script functions 

parameter_check $1 # either ID99 or ID95 for the 99/95% identity lastz data or Rm99 or Rm95 for the Repeatmasked LASTZ data at 99/95% identity 

create_folder tsv_files
create_folder blast.tables

for f in *.out.tsv #all files from perl extraction of tsv files containing FASTA DNA sequences from the pig X chromosome within the LASTZ hit start and stop positions
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

 for f in *.out.tsv.*; do mv "$f" "${f/.out.tsv./.}"; done #remove .out.tsv from file name (shortens the file names)

create_folder refined.blast
mv *.final refined.blast

mv *.out.tsv tsv_files
mv *.blast blast.tables

create_folder $1.blast
create_folder $1.blast/bar.graphs
create_folder $1.blast/split.files
create_folder $1.blast/subject.tables
create_folder $1.blast/subject.chromosome.positions
create_folder $1.blast/subject.tsv.files
create_folder $1.blast/subjects
create_folder $1.frequency.no.lcomplex

for f in refined.blast/$1**.blast.out.final; do cut -f8,9 $f >> $1.intermediate; done #for all refined blast files in either the ID99, Rm99, ID95 or Rm95 datasets extract the subject accession and title 
grep -v "chromosome" $1.intermediate > $1.intermediate2 # remove subject sequences containing the word chromosome as these sequences are low complexity regions detected in chromosome assemblies 
grep -v "scaffold" $1.intermediate2 > $1.titles.accessions # remove subject sequences containing the word scaffold as these sequences are low complexity regions detected in scaffold assemblies 
cut -f2 $1.titles.accessions > $1.accessions #extract the accession numbers from the data 
rm $1.intermediate
rm $1.intermediate2

/home/sq16564/Analysis-of-pig-X/./frequency.tables.R $1.accessions #Pass the accession numbers through the frequency.tables.R script to calculate the frequency of unique accession numbers in decending order 

head -n 20 $1.accessions.freq > top.20.acc.$1 #extract the top 20 most frequent accession numbers to isolate the most frequent subject sequences as most blast results have hundreds of returned sequences 

/home/sq16564/Analysis-of-pig-X/./frequency.bar.graph.R top.20.acc.$1 #create a frequency bar graph of the top 20 most frequent subject sequences 
create_folder $1.frequency.no.lcomplex/bar.graphs
mv *.png $1.frequency.no.lcomplex/bar.graphs

rm top.20.acc.$1

tail -n +20 $1.accessions.freq > remaining.acc.$1 #after extracting the top 20 most frequent subject sequences extract the remaining accessions 

split -l 50 remaining.acc.$1 #split the remaining accessions into smaller files as there are too many hits to plot into one graph
for f in x*; do log "creating frequency bar graph $f"; /home/sq16564/Analysis-of-pig-X/./frequency.bar.graph.R $f; done #for each split file create a frequency bar graph
rm remaining.acc.$1
create_folder $1.frequency.no.lcomplex/split.files
mv *.png $1.frequency.no.lcomplex/bar.graphs
mv x* $1.frequency.no.lcomplex/split.files

mv $1.accessions $1.frequency.no.lcomplex/
mv $1.titles.accessions $1.frequency.no.lcomplex/


awk -v freq=$2 '(NR==1) || ($2 > freq ) ' $1.accessions.freq > highest.freq.$1 # filter the subject sequences by the number of occurances to remove any sequences found once or twice where there are sequences with thousands of hits
# freq is $2 passed to the script as the identity and masking of a dataset will largely increase of decrease the number of occurances and the filtering must be edited accordingly

cut -f1 highest.freq.$1 > highest.accessions.$1 #extract the highest frequency accessions filtered by number of occurances from the datatframe 

rm highest.freq.$1

cat refined.blast/$1**.final >> $1.combined.blast #combine all the refined blast tables for the given dataset in $1 

cut -f8,9 $1.combined.blast > $1.combined.blast.titles #extract the subject title and accessions from the combined refined blast tables to have a list of subjects and accessions together 

rm $1.combined.blast

awk  '{ if (!a[$3]++ ) print ;}' $1.combined.blast.titles > $1.unique.blast.titles #create a list of the unique blast titles as the previous list of titles and accessions contains duplicates 

rm $1.combined.blast.titles

grep -f highest.accessions.$1 $1.unique.blast.titles > $1.most.frequent.blast.titles #create a list of the most frequent blast titles by extracting the highest accessions from the unique blast subject sequences file 


while read -r line #for every accession number in the highest frequency accessions
do
for f in refined.blast/$1**.final; do grep $line $f >> $line; done #take the accession and extract the blast table data from each blast results table  for the given accession 
cut -f1,4,5,6,7 $line > $line.position #extract the lastz hit information, the query start and stop, and the subject start and stop position 
sed 's/:/\t/g' $line.position > tab #separate the LASTZ hit data from : separation to tab separation 
cut -f5,6,8,9,10,11 tab > $line.tab #extract the lastz hit start and stop, the query start and stop, and the subject start and stop positions
/home/sq16564/Analysis-of-pig-X/./blast.subject.tables.R $line.tab #calculate the chromosome position of the subject sequence through the blast.subject.tables.R function
rm $line.position
rm tab
split -l 100 $line.tab.chromosome.hits $line.chromosome. #split the subject position table into smaller files 
for l in $line.chromosome.*; do perl /home/sq16564/Analysis-of-pig-X/./perl.extract.tsv.files.pl $l;done #extract the FASTA DNA sequence for the subject positions for manual blasting
mv $line $1.blast/subjects
mv $line.tab $1.blast/subject.tables
mv $line.tab.blast.hits $1.blast/subject.tables
mv $line.tab.chromosome.hits.out.tsv $1.blast/subject.tsv.files
mv $line.tab.chromosome.hits* $1.blast/subject.chromosome.positions
done < highest.accessions.$1


mv *$1* $1.blast
mv $1.accessions.freq $1.frequency.no.lcomplex/
