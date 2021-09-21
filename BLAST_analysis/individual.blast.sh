#!/bin/bash

#$ -cwd
#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas

#The following script is similar to full.blast.filtering.sh however it does not perform the initial blast search, it uses the blast data to produce output based on the subjects in the blast results table
#the input applies to individual BLASTED genes and their identity and masking properties 
source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #call the script containing the shell script functions 

parameter_check $1 #gene to be investigated and the identity and masking e.g ID99ENSSSCG000000 (99% identity LASTZ hits overlapping the gene name of with the ensembl ID)

create_folder $1.blast
create_folder $1.blast/bar.graphs
create_folder $1.blast/split.files
create_folder $1.blast/subject.tables
create_folder $1.blast/subject.chromosome.positions
create_folder $1.blast/subject.tsv.files
create_folder $1.blast/subjects
create_folder $1.blast/subject.blast

for f in refined.blast/$1**.blast.out.final; do cut -f8,9 $f >> $1.intermediate; done #for all refined blast files in either the ID99, Rm99, ID95 or Rm95 datasets extract the subject accession and title 
grep -v "chromosome" $1.intermediate > $1.intermediate2  # remove subject sequences containing the word chromosome as these sequences are low complexity regions detected in chromosome assemblies 
grep -v "scaffold" $1.intermediate2 > $1.titles.accessions # remove subject sequences containing the word scaffold as these sequences are low complexity regions detected in scaffold assemblies 
cut -f2 $1.titles.accessions > $1.accessions #extract the accession numbers from the data 
rm $1.intermediate
rm $1.intermediate2

/home/sq16564/Analysis-of-pig-X/./frequency.tables.R $1.accessions #Pass the accession numbers through the frequency.tables.R script to calculate the frequency of unique accession numbers in decending order 

split -l 50 $1.accessions.freq #split the accessions into smaller files as there are too many hits to plot into one graph
for f in x*; do log "creating frequency bar graph $f"; /home/sq16564/Analysis-of-pig-X/./frequency.bar.graph.R $f; done #for each split file create a frequency bar graph
rm remaining.acc.$1
create_folder $1.blast/split.files
mv *.png $1.blast/bar.graphs
mv x* $1.blast/split.files


awk '(NR==1) || ($2 > 1 ) ' $1.accessions.freq > highest.freq.$1 # filter the subject sequences by the number of occurances to remove any sequences found once or twice where there are sequences with thousands of hits
# freq is $2 passed to the script as the identity and masking of a dataset will largely increase of decrease the number of occurances and the filtering must be edited accordingly

cut -f1 highest.freq.$1 > highest.accessions.$1 #extract the highest frequency accessions filtered by number of occurances from the datatframe 

rm highest.freq.$1

cat refined.blast/$1**.final >> $1.combined.blast  #combine all the refined blast tables for the given dataset in $1

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
mv $line.tab.blast.hits $1.blast/subject.blast
mv *.tsv $1.blast/subject.tsv.files
mv $line.*chromosome.*  $1.blast/subject.chromosome.positions
done < highest.accessions.$1


mv *$1* $1.blast
mv $1.titles.accessions $1.blast
