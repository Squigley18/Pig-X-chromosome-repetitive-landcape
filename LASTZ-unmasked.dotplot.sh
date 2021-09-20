#!/bin/bash

#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas


source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh"
#pull functions from function script 


parameter_check $1 #check parameters have been passed 

touch lock.U.$2.r
 
lastz $1 --self --nomirror --notrivial --hspthresh=32500 --identity=$2 --format=rdotplot > $3.R #lastz command runs alignment on $1
#--self performs a self alignment
#--nomirror prevents mirror image alignments being shown
#--notrivial prevents output of trivial self alignment block 
#--hspthresh Sets score threshold for the x-drop extension method; HSPs scoring lower are discarded.
#--identity sets the percentage identity at which to filter the hits
#--format=dotplot creates data frame which can be used to create the naieve alignment dotplot


rm  lock.U.$2.r
