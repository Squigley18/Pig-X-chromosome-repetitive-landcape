#!/bin/Rscript
args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 


if (length(args)==0) { # test if there is at least one argument: if not, return an error
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) { #if there is a parameter passed to script then do the following

dots <- read.table(args[1],header=T) #read LASTZ dotplot table passed to script 

png(paste0(args[1],".png")) #save as png file with same name as LASTZ dotplot file .png
plot(dots, type="l") #create dot plot
dev.off()

}
