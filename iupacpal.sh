#!/bin/bash
# File: PlotWrapper.sh
#
#$ -cwd
#$ -S /bin/bash
#$ -m be
#$ -M sq16564@essex.ac.uk
#$ -q all.q

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #Call shell script functions from pre-defined script 


IUPACpal -f $1 -s $2 -m 50 -M 1000 -g 100 -x 0 -o $3.iupac
#To run the inverted repeat analysis the IUPACpal function must be called
#The DNA sequence is then defined using the -f parameter 
#The sequence header must also be provided to run the IUPACpal software using the -s parameter
#The minimum length of the palindromes to be detected is defined, in this case 50bp which removes the liklihood of detecting low complexity DNA using the -m parameter
#The minimum length of the palindromes to be detected is defined, in this case 1000bp allows for encompassing the average length of long palindromes found in other species but also reduced the computational demand using the -M parameter
#The maximum gap length was set as the default value of 100bp using the -g parameter
#The number of mismatched was set as the default value of 0bp using the -x parameter 
#Finally the text output is saved to the given output file defined with the -o parameter
