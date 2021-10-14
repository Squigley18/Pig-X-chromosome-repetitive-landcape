
#!/bin/bash

#$ -q all.q
#$ -M sq16564@essex.ac.uk
#$ -m beas


source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #pull functions from functions file


lastz $1[94000000..95000000]  --self --nomirror --notrivial --notransition --nogapped --chain --step=100  --identity=$2  --format=rdotplot > $3


