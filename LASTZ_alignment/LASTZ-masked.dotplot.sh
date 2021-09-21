#!/bin/bash

#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas


source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh"
#pull functions from function script 


parameter_check $1 #check parameters have been passed to script 

touch lock.M.$2.r
 
lastz $1 --self --nomirror --notrivial --identity=$2 --format=rdotplot > $3.R #lastz command runs alignment on $1
#--self performs a self alignment
#--nomirror prevents mirror image alignments being shown
#--notrivial prevents output of trivial self alignment block 
#--identity sets the percentage identity at which to filter the hits
#--format=dotplot creates data frame which can be used to create the naieve alignment dotplo

rm  lock.M.$2.r
