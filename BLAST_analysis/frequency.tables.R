#!/bin/Rscript

#librarys used in R analysis

library(tidyverse) #Tidyverse library is needed for the use of functions such as ggplot
library(GenomicFeatures) #GenomicFeatures library is required to load the TXdb object for the segment graphs 
library(GenomicRanges) #GenomicRanges library is required to create GenomicRanges objects 
library(ggbio) #Ggbio library is required for annotating the genes within the segment graph 
library(biomaRt) #Biomart library is required for accessing the gene data of the pig X chromosome 

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 


if (length(args)==0) { # test if there is at least one argument: if not, return an error
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) { #if there is a parameter passed to script then do the following
 
 # Read table of subject sequence accession numbers into R 
 query.subject <- read.table(args[1], header = FALSE, sep = "", 
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE)

#Calculate the frequency of occurances for each unique accession number to determine the most common subject sequences                            
frequency.table<- as.data.frame(sort(table(query.subject),decreasing=TRUE))

#write table of unique accession numbers and the number of occurances for creating a frequency bar graph and to allow for focus on the highest frequency results 
write.table(frequency.table,paste0(args[1],".freq"),sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
}

