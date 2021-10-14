#!/bin/bash
# File: PlotWrapper.sh
#
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

./LASTZ-dotplot.R $1

