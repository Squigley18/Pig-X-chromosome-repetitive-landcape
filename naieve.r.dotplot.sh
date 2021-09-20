#!/bin/bash
# File: PlotWrapper.sh
#
#$ 
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

source "/home/sq16564/Analysis-of-pig-X/shell-script-functions.sh" #access shell script functions within script shell-script-functions.sh


parameter_check $1

while [ -e $1 ] || [ -e $2 ] || [ -e $3 ] || [ -e $4 ]
 do
        log "sleeping until $1 $2 $3 $4 removed"

if [ -e $1 ]; then
  log "$1 still exists"
fi

if [ -e $2 ]; then
   log "$2 still exists"
fi

if [ -e $3 ]; then
  log "$3 still exists"
fi

if [ -e $4 ]; then
   log "$4 still exists"
fi

sleep 1800
done  #while loop waiting for lock files created in LASTZ-unmasked.dotplot.sh and LASTZ-masked.dotplot.sh to be removed before script continues 


create_folder naieve.r.alignment

/home/sq16564/Analysis-of-pig-X/LASTZ-analysis/./LASTZ-dotplot.R identity_99.R #create naieve dotplot for 99% identity LASTZ hits 
/home/sq16564/Analysis-of-pig-X/LASTZ-analysis/./LASTZ-dotplot.R identity_95.R #create naieve dotplot for 95% identity LASTZ hits 
/home/sq16564/Analysis-of-pig-X/LASTZ-analysis/./LASTZ-dotplot.R repeatmasked_99.R #create naieve dotplot for 99% identity repeatmasked LASTZ hits 
/home/sq16564/Analysis-of-pig-X/LASTZ-analysis/./LASTZ-dotplot.R repeatmasked_95.R #create naieve dotplot for 99% identity repeatmasked LASTZ hits 

mv *.R.png naieve.r.alignment
