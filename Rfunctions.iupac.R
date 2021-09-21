#!/bin/Rscript

#librarys used in R analysis
library(tidyverse) #Tidyverse library is needed for the use of functions such as ggplot
library(GenomicFeatures) #GenomicFeatures library is required to load the TXdb object for the segment graphs 
library(GenomicRanges) #GenomicRanges library is required to create GenomicRanges objects 
library(ggbio) #Ggbio library is required for annotating the genes within the segment graph 
library(biomaRt) #Biomart library is required for accessing the gene data of the pig X chromosome 

# Load the palindromes into an R readable table

# The load palindrome function returns a data frame containing the following columns:
# 'chromosome','y.col','start','stop','start2','end2',and 'width'

# palindrome - iupac palindrome file from bash


load.palindrome <- function(palindrome) {
  align <- read.table(palindrome, header = TRUE, sep = "",
                      numerals = c("allow.loss", "warn.loss","no.loss"),
                      col.names = c("start","stop","start2","stop2"),
                      stringsAsFactors = TRUE)

  data.table <- data.frame(align$start,align$stop,align$start2,align$stop2)
  data.table <- data.table %>% dplyr::rename(start = align.start,stop = align.stop,
                                             start2 = align.start2,stop2 = align.stop2)
  strand <- "*"
  chromosome <- "X"
  hitno <- cbind(chromosome,strand,data.table)
}

# GRfromDF function provides pre-defined arguments to the makeGRangesfromDataFrame function

# GRfromDF returns a Genomic Ranges dataframe containing the following columns:
# 'chromosome', 'IRanges', 'strand', and no metadata columns

# hit - the object name for the load.palindrome objects

GRfromDF <- function(hit){
makeGRangesFromDataFrame(hit,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("chromosome"),
                         start.field="start",
                         end.field=c("stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE )  
                             }
                             
#GRfromDF2 a similar function to GFfromDF but for creating Genomic Ranges using the start2 and end2 columns from the data this is to utilise the reciprical palindrome hits for the start and stop columns

# GRfromDF2 also  returns a Genomic Ranges dataframe containing the following columns: 
# 'chromosome', 'IRanges', 'strand', and no metadata columns 

# hit - the object name for the load.palindrome objects

GRfromDF2 <- function(hit){
  makeGRangesFromDataFrame(hit,
                           keep.extra.columns=FALSE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("chromosome"),
                           start.field="start2",
                           end.field=c("stop2"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE) 
}         

# Create Genomic Ranges dataframe of the chosen window size to investigate

# make.windows returns the Genomic Ranges dataframe ranging from 1 to the chromosome length (126e6 for the pig X chromosome) 

# window.size - desired size of each range
# chr.size - size of the chromosome investigated (126e6 for pig X chromosome)

make.windows = function(window.size, chr.size){
  starts = seq(1, chr.size-window.size+1, window.size)
  ends   = seq(window.size, chr.size, window.size)
  GRanges(seqnames="X",
          ranges=IRanges(
            start= starts,
            end  = ends),
          seqlengths=c(X=chr.size))
}

# The average lengths of the lastz hits led to the decision a window size of 3000 would cover the hits and for continuity the window size of 3000 was maintained, the window was created as follows: 

window.3000 <- make.windows(3000,126e6)

# The following function was created to take the input and return a dot plot with the X chromosome plot on the Xaxis and the number of palindrome hits on the Y axis

# df - the data frame made from the countOverlaps function

# mytitle - title to be given to the graph

overlap.hits.to.chromosome <- function(df,mytitle){
  overlap.hit.plot <- ggplot(df,aes(x=xaxis,na.rm=TRUE))+
    geom_point(aes(y=yaxis),na.rm=TRUE)+
    scale_x_continuous(name="Location on X chromosome (mbp)", label=c(0,30,60,90,120), breaks=c(0,10000,20000,30000,40000))+
    theme_classic()+
    ylab("Number of palindromes")+
    labs(title=mytitle)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

# Fetch genes on the X chromosome from the given regions.

# starts - an integer vector of start coordinates
# stops - an integer vector of stop coordinates

# fetch.x.genes.in.region returns a dataframe with the following columns:
# 'chromosome_name','ensembl_gene_id','external_gene_name',
# 'gene_biotype','start_position','end_position'

fetch.x.genes.in.region <- function(starts, stops){
  ensembl <- useMart("ensembl")
  ensembl    <- useDataset("sscrofa_gene_ensembl", mart = ensembl)
  filt       <- list(chromosome="X",start=starts,stop=stops)
  getBM(attributes=c('chromosome_name','ensembl_gene_id','external_gene_name',
                     'gene_biotype','start_position','end_position'),
                  filters = c('chromosome_name','start','end'),
                  values=filt,
                  mart = ensembl)
}

# Count overlaps between the palindrome hits and the genes within the selected region

# x.genes.in.region - the output from the fetch.x.genes.in.region function
# pal.GR - the genomic ranges object of the palindrome 
# The count.overlaps.in.region function returns a dataframe with the following columns:
#'region.overlap.pal','xaxis'

count.overlaps.in.region <- function(x.genes.in.region,pal.GR,m.pal.GR){
region.GR <- makeGRangesFromDataFrame(x.genes.in.region, keep.extra.columns=FALSE, ignore.strand=FALSE,
seqinfo=NULL, seqnames.field=c("chromosome_name"),
 start.field="start_position", end.field=c("end_position"), starts.in.df.are.0based=FALSE)

region.overlap.pal <- countOverlaps(region.GR,pal.GR)
region.overlap.pal <- as.data.frame(region.overlap.pal)

m.region.overlap.pal <- countOverlaps(region.GR,m.pal.GR)
m.region.overlap.pal <- as.data.frame(m.region.overlap.pal)

xaxis <- seq(1:nrow(region.overlap.pal))
region.overlaps <- cbind(xaxis,region.overlap.pal,m.region.overlap.pal)
}

# Create genomic ranges object for the gene start and stop position

# start_position - gene start position
# Stop_position - gene stop position

create.gene.gr <- function(start_position,stop_position){
  gene.gr <- GRanges("X", IRanges(start_position, stop_position))
}

# Create table with the gene start and stop position, the chromosome and value for y axis in graphs (y.column)

# start_position - gene start position

# Stop_position - gene stop position

create.gene.table <- function(start_position,stop_position){
  y.col <- seq(1)
  chromosome <- "X"
  start <- start_position
  stop <- stop_position
  gene <- data.frame(cbind(y.col,chromosome,start,stop))
  gene$y.col <- as.integer(as.character(gene$y.col))
  gene$chromosome <- as.character(gene$chromosome)
  gene$start <- as.integer(as.character(gene$start))
  gene$stop <- as.integer(as.character(gene$stop))
  gene
}

# Create a TxDb object for the sus scrofa dataset (required for some of the following analysis)

Txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                    dataset="sscrofa_gene_ensembl", transcript_ids=NULL,
                    filter=NULL,
                    id_prefix="ensembl_", host="www.ensembl.org",
                    port=80, taxonomyId=NA, miRBaseBuild=NA)

# Create table containing the chromosome, width, start and stop positions of the palindrome hits matching within the gene of interest

# start_position - gene start position
# Stop_position - gene stop position
# mydf - the palindrome table data frame

subset.overlaps.for.gene.of.interest <- function(mydf,start_position,stop_position){
  gene.gr <- create.gene.gr(start_position,stop_position)
  gene.table <- create.gene.table(start_position,stop_position)
  my.data.frame.gr <- GRfromDF(mydf)
  find.overlaps <- findOverlaps(query = gene.gr, subject = my.data.frame.gr, type = 'any')
  overlap.df <-  data.frame(mydf[subjectHits(find.overlaps),], gene.table[queryHits(find.overlaps),])
  final.data.set <- overlap.df[,c(1,2,3,4,7)]
}

# Create table containing the chromosome, width, start and stop positions of the palindrome hits matching within the gene of interest

# start_position - gene start position

# Stop_position - gene stop position

# mydf - the palindrome table data frame

subset.overlaps.for.tables <- function(mydf,start_position,stop_position){
    gene.gr <- create.gene.gr(start_position,stop_position)
    gene.table <- create.gene.table(start_position,stop_position)
    my.data.frame.gr <- GRfromDF(mydf)
    find.overlaps <- findOverlaps(query = gene.gr, subject = my.data.frame.gr, type = 'any')
    overlap.df <-  data.frame(mydf[subjectHits(find.overlaps),], gene.table[queryHits(find.overlaps),])
    final.data.set <- overlap.df[,c(1,2,3,6)]
}

# Using the output from the subset.overlaps.for.tables function, the subset.overlaps.for.graphs function adds a column for the y axis (y.column) needed for some of the graphs

# start_position - gene start position

# Stop_position - gene stop position

# mydf - the palindrome table data frame

subset.overlaps.for.graphs <- function(mydf,start_position,stop_position){
  subset.table <- subset.overlaps.for.gene.of.interest(mydf,start_position,stop_position)
  for.graph.data <- cbind(subset.table,y.column = seq(from= 2,by= 1,length.out=nrow(subset.table)))
}

# Plot the palindromes as segments overlapping the ggbio annotated gene

# start_position - gene start position

# Stop_position - gene stop position

# mydf - dataframe created using the subset.overlaps.for.graphs function

# mytitle - a character string representing the title of the plot

# subset.overlap.segment.plot returns a track object which can be used to view the ggplot created or save the ggplot using the ggsave function

subset.overlap.segment.plot <- function(mydf,start_position,stop_position,mytitle) {
  gene.gr <- create.gene.gr(start_position,stop_position)
  gene <- create.gene.table(start_position,stop_position)
  subset.data <- subset.overlaps.for.graphs(mydf,start_position,stop_position)
  overlap.plot <<- autoplot(Txdb, which = gene.gr, label.color="grey20") +
    theme_classic()+
    ylab("Palindromes")+
    geom_segment(data=subset.data,aes(x=start,xend=stop,y=y.column,yend=y.column),color='blue')+
    geom_segment(data=gene,aes(x=start,xend=stop,y=y.col,yend=y.col),color='red')+
    theme(axis.title.x = element_text(color = "grey20", size = 9, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 8, angle = 0, hjust = .5, vjust = 0, face = "bold"))
  track <<- tracks(overlap.plot,title=mytitle,xlab="Location on the chromosome")}
  
# Create table containing the chromosome, width, start2 and stop2 positions of the  palindrome hits within the gene of interest 

# position 2 is referring to the use of start2 and stop2 from the table compared to start1 and stop1 

# start_position - gene start position 

# Stop_position - gene stop position 

# mydf - the palindrome table data frame 

subset.overlaps.for.gene.of.interest.position.2 <- function(mydf,start_position,stop_position){
  gene.gr <- create.gene.gr(start_position,stop_position)
  gene.table <- create.gene.table(start_position,stop_position)
  my.data.frame.gr <- GRfromDF(mydf)
  find.overlaps <- findOverlaps(query = gene.gr, subject = my.data.frame.gr, type = 'any')
  overlap.df <-  data.frame(mydf[subjectHits(find.overlaps),], gene.table[queryHits(find.overlaps),])
}

#Using the subset.overlaps.for.gene.of.interest.position.2 function the subset hits 
#can be used to create a data frame to overlap to the ideogram to show the distribution of palindrome hits 
#Where there are 0 overlaps detected the values are converted to NAs in the ideogram column and the hits above 1 are converted to the corresponding xaxis row value

#mydf- the output from subset.overlaps.for.gene.of.interest.position.2

manipulate.data.for.ideogram.plot <- function(mydf){
  hit.gr <- GRfromDF2(mydf)
  window.3000 <- make.windows(3000,126e6)
  hits <- countOverlaps(window.3000,hit.gr)
  xaxis <- seq(1:42000)
  counts <- as.data.frame(cbind(xaxis,hits))
  counts$ideogram <- counts$hits
  counts$ideogram[counts$hits >= 1] <- counts$xaxis[counts$hits >= 1]
  counts$ideogram[counts$ideogram==0] <- NA
  counts
}


# The following function is used to produce an ideogram plot with two ideograms of the pig X chromosome where the repeat masked hits are on the top ideogram and the unmasked hits on the bottom ideogram.

# df- is the data frame created using the above function manipulate.data.for.ideogram

# mytitle - the title given to the ideogram plot

Ideogram.plot <- function(mydf,mytitle){
  chromosome.ideogram <- data.frame(x = 1:42000, y = c(rep.int(1,20500),
                                                       rep.int(NA,800),rep.int(1,20700)),
                                                     z = c(rep.int(2,20500),
                                          rep.int(NA,800),rep.int(2,20700)))
  
  Ideogram <- ggplot(chromosome.ideogram, aes(x,y),na.rm = TRUE)+ 
    geom_path(size = 4, lineend = "round",colour="gray87")+
    geom_path(aes(x,z),size = 4, lineend = "round",colour="gray87")+
    geom_segment(mydf,aes(x=palindrome,xend=palindrome+15,
                                   y=1,yend=1),size=3,colour="black")+
    geom_segment(mydf,aes(x=m.palindrome,xend=m.palindrome+15,
                                   y=2,yend=2),size=3,colour="black")+
    scale_x_continuous(name="Location on chromosome (mbp)", label=c(0,30,60,90,120), breaks=c(0,10000,20000,30000,40000))+
    scale_y_discrete(na.omit(chromosome.ideogram,FALSE))+
    annotate("text", x = c(1,1.5), y = c(1.5,2.5),
             label = c("Unmasked", "Repeatmasked") , color="black",
             size=1.5 , fontface="bold")+
    theme_classic()+
    labs(title=mytitle)+
    theme(axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),aspect.ratio = 1/8)
}
