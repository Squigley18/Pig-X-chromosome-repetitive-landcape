#!/bin/bash 

#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas


source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #pull functions from functions file

parameter_check $1 #function to check parameters have been passed to script, if not the script is killed
file_check $3 #function to check the file $3 does not exit otherwise script is killed 

create_folder alignment_$4
touch lock.M.$2

lastz $1 --self --nomirror --notrivial  --identity=$2 --format=general:strand1,start1,\
end1,strand2,start2,end2,identity > $3 #lastz command runs alignment on $1
#--self performs a self alignment
#--nomirror prevents mirror image alignments being shown
#--notrivial prevents output of trivial self alignment block 
#--identity sets the percentage identity at which to filter the hits 
#--format flag with chosen parameters provides a table showing the; start position, stop position and strand for the query sequence;
# start position, stop position and strand for the target sequence; the percentage identity between the sequences, and the coverage between the sequences 

mv $3 alignment_$4
rm  lock.M.$2
