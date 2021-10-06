#!/bin/Rscript
library(tidyverse) #Tidyverse library is needed for the use of functions such as ggplot
library(GenomicFeatures) #GenomicFeatures library is required to load the TXdb object for the segment graphs 
library(GenomicRanges) #GenomicRanges library is required to create GenomicRanges objects 
library(ggbio) #Ggbio library is required for annotating the genes within the segment graph 
library(biomaRt) #Biomart library is required for accessing the gene data of the pig X chromosome 

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 


if (length(args)==0) { # test if there is at least one argument: if not, return an error
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) { #if there is a parameter passed to script then do the following
 
# Read table of subject sequence accession numbers and their frequency count from frequency.tables.R script output 
sub <- read.table(args[1], header = FALSE, sep = "", 
                      numerals = c("allow.loss", "warn.loss","no.loss"),
                      stringsAsFactors = TRUE) 
 
# plot the frequency counts for each accession number as a bar graph for comparison 
plot <- ggplot(sub, aes(x=reorder(V1,-V2),y=V2)) +
  geom_col()+
  ylab("Number of occurences")+
  xlab("Accession number")+
  labs(title="frequency of  BLAST results")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,size=5),
        axis.title.x = element_text(size = 10),
		axis.title.y = element_text(size = 10))
  ggsave(paste0(args[1],".bar.svg"),plot,unit="mm",dpi=300,height=85,width=120)

}
