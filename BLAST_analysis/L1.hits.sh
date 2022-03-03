#!/bin/bash

#$ -cwd
#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas

#The following script is similar to full.blast.filtering.sh however it does not perform the initial blast search, it uses the blast data to produce output based on the subjects in the blast results table
#the input applies to individual BLASTED genes and their identity and masking properties
source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #call the script containing the shell script functions

grep "L1" refined.blast/Rm99**.final > Rm99.L1
grep "L1" refined.blast/Rm95**.final > Rm95.L1

cut -f1,4,5,6,7 Rm99.L1 > Rm99.positions
cut -f1,4,5,6,7 Rm95.L1 > Rm95.positions

sed 's/:/\t/g' Rm99.positions > Rm99.tab
sed 's/:/\t/g' Rm95.positions > Rm95.tab

cut -f6,7,9,10,11,12 Rm99.tab > Rm99.tab2
cut -f6,7,9,10,11,12 Rm95.tab > Rm95.tab2

/home/sq16564/Analysis-of-pig-X/./blast.subject.tables.R Rm99.tab2
/home/sq16564/Analysis-of-pig-X/./blast.subject.tables.R Rm95.tab2

/home/sq16564/./L1.ideo.R Rm99.tab2.chromosome.hits Rm95.tab2.chromosome.hits
