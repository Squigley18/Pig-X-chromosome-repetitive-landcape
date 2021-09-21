#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 
  
source("/home/sq16564/Analysis-of-pig-X/Rfunctions.R")
  
ID99 <- load.alignment(args[1],99) #load 99% identity LASTZ hit table into R 
ID95 <- load.alignment(args[2],99) #load 95% identity LASTZ hit table into R 
IDrm99 <- load.alignment(args[3],99) #load 99% identity repeatmasked LASTZ hit table into R 
IDrm95 <- load.alignment(args[4],99) #load 95% identity repeatmasked LASTZ hit table into R 

#Create a dataframe of the combined alignments which will be used when looping certain actions later

a <- data.frame(group = "ID95", value = ID95)
b <- data.frame(group = "ID99", value = ID99)
c <- data.frame(group = "Rm99", value = IDrm99)
d <- data.frame(group = "Rm95", value = IDrm95)
combined.alignments <- rbind(a,b,c,d)
combined.alignments <- combined.alignments %>% dplyr::rename(chromosome=value.chromosome,
                                                             start = value.start, stop=value.stop,
                                                             start2 = value.start2,
                                                             end2 = value.end2,width=value.width)

  
ID99gr <- GRfromDF(ID99) #create genomic ranges object for 99% identity LASTZ hits start and stop positions 
ID95gr <- GRfromDF(ID95) #create genomic ranges object for 95% identity LASTZ hits start and stop positions 
IDrm99gr <- GRfromDF(IDrm99) #create genomic ranges object for 99% identity repeatmasked LASTZ hits start and stop positions 
IDrm95gr <- GRfromDF(IDrm95) #create genomic ranges object for 95% identity repeatmasked LASTZ hits start and stop positions 
  
#load the Txdb object from the functions script
Txdb

# fetch.x.genes.in.region returns a dataframe with the following columns:
# 'chromosome_name','ensembl_gene_id','external_gene_name',
# 'gene_biotype','start_position','end_position'

# genes are fetched from positions along length of the chromosome as follows excluding PAR (not of interest in this study at 0-7000000bp):

region1 <- fetch.x.genes.in.region ("7000000","27000000")
region2 <- fetch.x.genes.in.region ("27000001","47000000")
region3 <- fetch.x.genes.in.region ("47000001","67000000")
region4 <- fetch.x.genes.in.region ("67000001","87000000")
region5 <- fetch.x.genes.in.region ("87000001","107000000")
region6 <- fetch.x.genes.in.region ("107000001","127000000")

#The count.overlaps,in.region function between the genes within the chromosome region and the LASTZ hits is performed to determine the number of overlapping hits 

R1 <- count.overlaps.in.region(region1,IDrm99gr,IDrm95gr,ID99gr,ID95gr)
R2 <- count.overlaps.in.region(region2,IDrm99gr,IDrm95gr,ID99gr,ID95gr)
R3 <- count.overlaps.in.region(region3,IDrm99gr,IDrm95gr,ID99gr,ID95gr)
R4 <- count.overlaps.in.region(region4,IDrm99gr,IDrm95gr,ID99gr,ID95gr)
R5 <- count.overlaps.in.region(region5,IDrm99gr,IDrm95gr,ID99gr,ID95gr)
R6 <- count.overlaps.in.region(region6,IDrm99gr,IDrm95gr,ID99gr,ID95gr)

#The next stage is to filter the genes so only the genes with over 0 hits overlapping the 99% identity hits or 99% identity repeat masked hits remain

r1.overlaps.99 <- R1 %>% filter(region.overlap.rm99 > 0 | region.overlap.99 > 0)
r2.overlaps.99 <- R2 %>% filter(region.overlap.rm99 > 0 | region.overlap.99 > 0)
r3.overlaps.99 <- R3 %>% filter(region.overlap.rm99 > 0 | region.overlap.99 > 0)
r4.overlaps.99 <- R4 %>% filter(region.overlap.rm99 > 0 | region.overlap.99 > 0)
r5.overlaps.99 <- R5 %>% filter(region.overlap.rm99 > 0 | region.overlap.99 > 0)
r6.overlaps.99 <- R6 %>% filter(region.overlap.rm99 > 0 | region.overlap.99 > 0)

#The values from the x axis column are then extracted as a list which is used to extract the filtered genes from the original region tables

r1.rows <- r1.overlaps.99$xaxis
r2.rows <-  r2.overlaps.99$xaxis
r3.rows <-  r3.overlaps.99$xaxis
r4.rows <-  r4.overlaps.99$xaxis
r5.rows <-  r5.overlaps.99$xaxis
r6.rows <-  r6.overlaps.99$xaxis

#The genes which have over 0 hits from the 99% identity and 99% identity repeat masked LASTZ tables overlapping them are extracted and held within these objects:

r1.genes <- region1[r1.rows,]
r2.genes <- region2[r2.rows,]
r3.genes <- region3[r3.rows,]
r4.genes <- region4[r4.rows,]
r5.genes <- region5[r5.rows,]
r6.genes <- region6[r6.rows,]

#All the extracted genes are combined into a table showing the following columns:
# 'chromosome_name','ensembl_gene_id','external_gene_name',
# 'gene_biotype','start_position','end_position'

all.gene <- rbind(r1.genes,r2.genes,r3.genes,r4.genes,r5.genes,r6.genes)

#The following loop takes the LASTZ hits at 99% identity and 95% identity both masked and unmasked
#For each gene in the all.gene table, 
#it takes the overlapping LASTZ hits and plots them overlapping the gene (annotated with ggbio) to determine the region of overlap and whether the overlap is intronic or exonix 
#The graphs are then saved using ggsave


  for(row in unique(all.gene$ensembl_gene_id)){
    sub <- all.gene[all.gene$ensembl_gene_id == row,c(5,6)]
    start <- sub$start_position
    stop <- sub$end_position
  for(group.name in unique(combined.alignments$group)){
  subset.data <- combined.alignments[combined.alignments$group == group.name,]
  sub.plot <- subset.overlap.segment.plot(subset.data,start,stop,paste(group.name,row))
  ggsave(paste0(group.name,row,"_segments.svg"),sub.plot,unit="mm",height=85,width=85)}}

#The following loop takes the LASTZ hits at 99% identity and 95% identity both masked and unmasked
#For each gene in the all.gene table, 
#it takes the overlapping LASTZ hits and plots their original locations on an ideogram to determine their distribution pattern throughout the chromosome
#The graphs are then saved using ggsave

for(rows in unique(all.gene$ensembl_gene_id)){ 
  sub <- all.gene[all.gene$ensembl_gene_id == rows,c(5,6)]
  start <- sub$start_position
  stop <- sub$end_position
  I99 <- subset.overlaps.for.gene.of.interest.position.2(ID99,start,stop)
  I95 <- subset.overlaps.for.gene.of.interest.position.2(ID95,start,stop)
  R99 <- subset.overlaps.for.gene.of.interest.position.2(IDrm99,start,stop)
  R95 <- subset.overlaps.for.gene.of.interest.position.2(IDrm95,start,stop)
  if (nrow(I99) == 0) { 
chromosome <- "X"
start2 <- 0
end2 <- 0
width <- 0
I99 <- as.data.frame(cbind(chromosome,start2,end2,width))
I99 }
  
  if (nrow(I95) == 0) { chromosome <- "X"
start2 <- 0
end2 <- 0
width <- 0
I95 <- as.data.frame(cbind(chromosome,start2,end2,width))
I95   }
  if (nrow(R99) == 0) { chromosome <- "X"
start2 <- 0
end2 <- 0
width <- 0
R99 <- as.data.frame(cbind(chromosome,start2,end2,width)) 
R99 }
  if (nrow(R95) == 0) { chromosome <- "X"
start2 <- 0
end2 <- 0
width <- 0
R95 <- as.data.frame(cbind(chromosome,start2,end2,width))
R95 }
  ID99.hits <- manipulate.data.for.ideogram.plot(I99)
  ID95.hits <- manipulate.data.for.ideogram.plot(I95)
  Rm99.hits <- manipulate.data.for.ideogram.plot(R99)
  Rm95.hits <- manipulate.data.for.ideogram.plot(R95)
  ideogram.hits <- as.data.frame(cbind(ID99.hits$xaxis,ID99.hits$ideogram,ID95.hits$ideogram,Rm99.hits$ideogram,Rm95.hits$ideogram))
  ideogram.hits <- ideogram.hits %>% dplyr::rename(Xaxis = V1, ID99 = V2, ID95 = V3, Rm99 = V4, Rm95 = V5)
  plot <- Ideogram.plot(ideogram.hits,paste("lastz hits showing homology to",rows,"original distribition"))
  ggsave(paste0(rows,"_ideogram.png"),plot,unit="mm",height=85,width=120)
}

#The following loop takes the LASTZ hits at 99% identity and 95% identity both masked and unmasked
#For each gene in the all.gene table, 
#it takes the overlapping LASTZ hits and exports the chromosome, start and stop positions in a table for extracting the DNA sequence using the perl script 


for(group.name in unique(combined.alignments$group)){
  subset.data <- combined.alignments[combined.alignments$group == group.name,]
  for(rows in unique(all.gene$ensembl_gene_id)){
    sub <- all.gene[all.gene$ensembl_gene_id == rows,c(5,6)]
    start <- sub$start_position
    stop <- sub$end_position
    table <- subset.overlaps.for.gene.of.interest(subset.data,start,stop)
    if( nrow(table) == 0){
        next
        }
    table.final <- table[,c(2,3,4)]
    names(table.final) <- NULL
    write.table(table.final,paste0(group.name,rows,".tab"),sep="\t",
                row.names = FALSE,quote = FALSE)}}
