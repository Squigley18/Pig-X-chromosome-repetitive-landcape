#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 

source('./Rfunctions.iupac.R') #Call the functions pre-defined in the script 
  
palindrome <- load.palindrome(args[1]) #Load the unmasked palindrome table of start and stop positions
m.palindrome <- load.palindrome(args[2]) #Load the masked palindrome table of start and stop positions

#Create a genomic ranges object for the palindromes

palgr <- GRfromDF(palindrome)
m.palgr <- GRfromDF(m.palindrome)
 
#Count the overlaps between the palindrome genomic ranges objects and the chromosome windows
overlap.pal <- countOverlaps(window.3000,palgr)
m.overlap.pal <- countOverlaps(window.3000,m.palgr) 

# The X chromosome GenomicRanges object has 42000 rows which represents the windows of 3000bp along the X chromosome 

xaxis <- seq(1:42000)

# Convert the countOverlaps objects into data frames which can then be used to create the dotplot through ggplot 

overlap.pal.df <- as.data.frame(cbind(xaxis,overlap.pal))
overlap.pal.df[overlap.pal.df==0] <- NA
overlap.pal.df <- overlap.pal.df %>% dplyr::rename(yaxis = overlap.pal)

m.overlap.pal.df <- as.data.frame(cbind(xaxis,m.overlap.pal))
m.overlap.pal.df[m.overlap.pal.df==0] <- NA
m.overlap.pal.df <- m.overlap.pal.df %>% dplyr::rename(yaxis = m.overlap.pal)


# The following are the objects for the dot plots created with the overlap.hits.to.chromosome function

pal.count.plot <- overlap.hits.to.chromosome(overlap.pal.df,"palindrome distribution")
m.pal.count.plot <- overlap.hits.to.chromosome(m.overlap.pal.df,"palindrome distribution")

ggsave(paste0(args[1],".overlap.svg"),pal.count.plot,dpi=300,unit="mm",height=85,width=120)
ggsave(paste0(args[2],".overlap.svg"),m.pal.count.plot,dpi=300,unit="mm",height=85,width=120)
 
#The following functions create the dataframe which can be passed to the ideogram plot function 
palindrome.hits <- manipulate.data.for.overlap.ideogram.plot(palindrome)
  m.palindrome.hits <- manipulate.data.for.overlap.ideogram.plot(m.palindrome)
#The dataframes are combined into one object to be passed to the ideogram plot function 
  all.ideogram.hits <- as.data.frame(cbind(palindrome.hits$xaxis,palindrome.hits$ideogram,m.palindrome.hits$ideogram))
  all.ideogram.hits <- all.ideogram.hits %>% dplyr::rename(Xaxis = V1, palindrome = V2, m.palindrome = V3)
#The ideogram plot can be created 
  all.plot <- Ideogram.overlap.plot(all.ideogram.hits,"all palindrome hits distribution")
  ggsave("all.palindromes.ideogram.svg",all.plot,unit="mm",height=85,width=120) #Saved ideogram plot

