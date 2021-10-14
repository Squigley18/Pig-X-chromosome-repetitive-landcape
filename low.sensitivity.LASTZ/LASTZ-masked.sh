
#!/bin/bash

#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas


source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #pull functions from functions file


lastz $1[94000000..95000000]  --self --nomirror --notrivial --notransition --nogapped --chain --step=100  --identity=$2  --format=rdotplot > $3
#lastz command runs alignment on $1 given regions
#--self performs a self alignment
#--nomirror prevents mirror image alignments being shown
#--notrivial prevents output of trivial self alignment block 
#--notransition reduces the sensitivity of the alignment
#--nogapped skips the gapped alignment stage of the LASTZ process
#--chain chains together high scoring segment pairs to form a high scoring segment path
#--step=100 reduces the sensitivity of the alignment
#--identity sets the percentage identity at which to filter the hits
#--format=dotplot creates data frame which can be used to create the naieve alignment dotplot
