#!/bin/bash
# File: PlotWrapper.sh
#
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #Call shell script functions pre-defined in script

parameter_check $1  #check $1 has been passed to script

#The DNA sequence was split in half to be passed through the iupacpal software separatly to reduce computational demand, the following code is used for both files for the unmasked DNA sequence

tail -n +14 $1 > iupac.out #remove header from iupacpal output text file (masked DNA iupacpal output file 1)
tail -n +14 $2 > m.iupac.out #remove header from iupacpal output text file (masked DNA iupacpal output file 2)

#Remove all lines from text file where N's are present due to masking or incomplete assemblies
awk 'BEGIN { RS="\n\n" } { if($0!~/nn/){print $0"\n"} }' iupac.out > refined.iupac 
awk 'BEGIN { RS="\n\n" } { if($0!~/nn/){print $0"\n"} }' m.iupac.out > refined.m.iupac

#Take every other line starting from line 2 which removes the alignment bars from the text file
sed -n '2~2!p' refined.iupac > condensed.iupac
sed -n '2~2!p' refined.m.iupac > condensed.m.iupac

rm iupac.out
rm m.iupac.out
rm refined.iupac
rm refined.m.iupac

#remove the white space and alignment sequences from the text file and extract the first column as the start position and the second column as the stop position 
cat condensed.iupac | tr '[a-z]' ' ' | tr -s ' ' | cut -d ' ' -f1 > start
cat condensed.m.iupac | tr '[a-z]' ' ' | tr -s ' ' | cut -d ' ' -f1 > m.start

cat condensed.iupac | tr '[a-z]' ' ' | tr -s ' ' | cut -d ' ' -f2 > stop
cat condensed.m.iupac | tr '[a-z]' ' ' | tr -s ' ' | cut -d ' ' -f2 > m.stop

rm condensed.iupac
rm condensed.m.iupac

#Take the alternate lines for the start and stop positions for the complimentary palindromes positions. 
sed -n '2~2!p' start > start1
sed -n '2~2!p' m.start > m.start1

sed -n '1~2!p' start > stop2
sed -n '1~2!p' m.start > m.stop2

sed -n '2~2!p' stop > stop1
sed -n '2~2!p' m.stop > m.stop1

sed -n '1~2!p' stop > start2
sed -n '1~2!p' m.stop > m.start2

rm start
rm stop 
rm m.start
rm m.stop

#Combine the seperate start and stop positions into one table with four columns; start1, stop1, start2, stop2
paste start1 stop1 start2 stop2 > $1.tab
paste m.start1 m.stop1 m.start2 m.stop2 > $2.tab

awk NF $1.tab #remove trailing white space
awk NF $2.tab #remove trailing white space

rm start1
rm start2
rm stop1
rm stop2
rm m.start1
rm m.start2
rm m.stop1
rm m.stop2

#Add headers to the columns to be used in a consistent format in the later analysis stages and combine the two unmasked DNA sequence files
./iupac.tables.R $1.tab $2.tab
