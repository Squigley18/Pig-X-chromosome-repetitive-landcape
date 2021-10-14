
#!/bin/bash

#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas


source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #pull functions from functions file


lastz $1[0..10000000] --self --nomirror --notrivial --notransition --nogapped --step=100 --hspthresh=31250 --identity=$2 --format=rdotplot > $3


