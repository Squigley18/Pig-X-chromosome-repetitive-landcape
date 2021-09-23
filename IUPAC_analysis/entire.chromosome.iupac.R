#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 
  
source("Rfunctions.iupac.R") #Call functions created in Rfunctions.iupac.R script
  
palindrome <- load.palindrome(args[1]) #load unmasked palindrome position table into R 
m.palindrome <- load.palindrome(args[2]) #load masked palindrome position table into R 


#Create a genomic ranges object for the palindromes

palgr <- GRfromDF(palindrome) #create genomic ranges object for unmasked palindrome hits start and stop positions 
m.palgr <- GRfromDF(m.palindrome) #create genomic ranges object for masked palindrome hits start and stop positions 
  
#load the Txdb object from the functions script
Txdb

# fetch.x.genes.in.region returns a dataframe with the following columns:
# 'chromosome_name','ensembl_gene_id','external_gene_name',
# 'gene_biotype','start_position','end_position'

# genes are fetched from length of the chromosome as follows excluding PAR (not of interest in this study at 0-7000000bp):

region1 <- fetch.x.genes.in.region ("7000000","27000000")
region2 <- fetch.x.genes.in.region ("27000001","47000000")
region3 <- fetch.x.genes.in.region ("47000001","67000000")
region4 <- fetch.x.genes.in.region ("67000001","87000000")
region5 <- fetch.x.genes.in.region ("87000001","107000000")
region6 <- fetch.x.genes.in.region ("107000001","127000000")

#The overlap count between the genes within the chromosome regions and the palindrome hits is performed

R1 <- count.overlaps.in.region(region1,palgr,m.palgr)
R2 <- count.overlaps.in.region(region2,palgr,m.palgr)
R3 <- count.overlaps.in.region(region3,palgr,m.palgr)
R4 <- count.overlaps.in.region(region4,palgr,m.palgr)
R5 <- count.overlaps.in.region(region5,palgr,m.palgr)
R6 <- count.overlaps.in.region(region6,palgr,m.palgr)

#The next stage is to filter the genes so only the genes with overlapping palindrome hits (masked and unmasked) remain

r1.overlaps <- R1 %>% filter(region.overlap.pal > 0 | m.region.overlap.pal > 0 )
r2.overlaps <- R2 %>% filter(region.overlap.pal > 0 | m.region.overlap.pal > 0 )
r3.overlaps <- R3 %>% filter(region.overlap.pal > 0 | m.region.overlap.pal > 0 )
r4.overlaps <- R4 %>% filter(region.overlap.pal > 0 | m.region.overlap.pal > 0 )
r5.overlaps <- R5 %>% filter(region.overlap.pal > 0 | m.region.overlap.pal > 0 )
r6.overlaps <- R6 %>% filter(region.overlap.pal > 0 | m.region.overlap.pal > 0 )

#The values from the x axis column are then extracted as a list which is used to extract the filtered genes from the original region tables

r1.rows <- r1.overlaps$xaxis
r2.rows <-  r2.overlaps$xaxis
r3.rows <-  r3.overlaps$xaxis
r4.rows <-  r4.overlaps$xaxis
r5.rows <-  r5.overlaps$xaxis
r6.rows <-  r6.overlaps$xaxis

#The genes which have over 0 hits palindromes overlapping them are extracted and held within these objects:

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

#The following loop takes the palindrome positions and for each gene in the all.gene.table
#it takes the overlapping palindromes and plots them overlapping the gene (annotated with ggbio) to determine where the overlap occurs and whether the overlap covers an intron or an exon 
#The graphs are then saved using ggsave

 for(row in unique(all.gene$ensembl_gene_id)){
    sub <- all.gene[all.gene$ensembl_gene_id == row,c(5,6)]
    start <- sub$start_position
    stop <- sub$end_position
  sub.plot <- subset.overlap.palindrome.plot(palindrome,start,stop,paste("palindromes aligned to", row))
  sub.plot.2 <- subset.overlap.palindrome.plot(m.palindrome,start,stop,paste("masked palindromes aligned to",row))
  ggsave(paste0(row,".pal_segments.svg"),sub.plot,unit="mm",height=85,width=85)
  ggsave(paste0(row,".m.pal_segments.svg"),sub.plot.2,unit="mm",height=85,width=85)}


#The following loop takes the palindromes within each gene in the all.gene table, 
#it takes the complimentary overlapping palindromes and plots their locations on an ideogram to determine the palindrome didstribution
#The graphs are then saved using ggsave

for(rows in unique(all.gene$ensembl_gene_id)){ 
  sub <- all.gene[all.gene$ensembl_gene_id == rows,c(5,6)]
  start <- sub$start_position
  stop <- sub$end_position
  pal <- subset.overlaps.for.gene.of.interest.position.2(palindrome,start,stop)
  mpal <- subset.overlaps.for.gene.of.interest.position.2(m.palindrome,start,stop)
  pal.hits <- manipulate.data.for.ideogram.plot(pal)
  mpal.hits <- manipulate.data.for.ideogram.plot(mpal)
  ideogram.hits <- as.data.frame(cbind(pal.hits$xaxis,pal.hits$ideogram,mpal.hits$ideogram))
  ideogram.hits <- ideogram.hits %>% dplyr::rename(Xaxis = V1, palindrome = V2, m.palindrome = V3)
  plot <- Ideogram.plot(ideogram.hits,paste("palindromes within",rows,"distribition"))
  ggsave(paste0(rows,".pal_ideogram.svg"),plot,unit="mm",height=85,width=120)
}

#The following loop takes the palindrome hits both masked and unmasked for each gene in the all.gene table, 
#it takes the overlapping palindrome hits and exports the chromosome, start and stop positions in a table


  for(rows in unique(all.gene$ensembl_gene_id)){ 
    sub <- all.gene[all.gene$ensembl_gene_id == rows,c(5,6)]
    start <- sub$start_position
    stop <- sub$end_position
    table <- subset.overlaps.for.gene.of.interest(palindrome,start,stop)
    table2 <- subset.overlaps.for.gene.of.interest(m.palindrome,start,stop)
    table.final <- table[,c(1,3,4)]
    table.final2 <- table2[,c(1,3,4)]
    names(table.final) <- NULL
    names(table.final2) <- NULL
    write.table(table.final,paste0(rows,".pal.tab"),sep="\t",
                row.names = FALSE,quote = FALSE)
    write.table(table.final2,paste0(rows,".m.pal.tab"),sep="\t",
                row.names = FALSE,quote = FALSE)
}


#Finally many of the palindromes had also not overlapped any genes in the X chromosome and therefore their complimentary palindromes were also plot to ideograms to determine their distribution
#These palindrome start and stop positions were exported into a table

  palindrome.hits <- manipulate.data.for.ideogram.plot(palindrome)
  m.palindrome.hits <- manipulate.data.for.ideogram.plot(m.palindrome)
  all.ideogram.hits <- as.data.frame(cbind(palindrome.hits$xaxis,palindrome.hits$ideogram,m.palindrome.hits$ideogram))
  all.ideogram.hits <- all.ideogram.hits %>% dplyr::rename(Xaxis = V1, palindrome = V2, m.palindrome = V3)
  all.plot <- Ideogram.plot(all.ideogram.hits,"all palindrome hits distribution")
  ggsave("all.palindromes.ideogram.svg",all.plot,unit="mm",height=85,width=120)

 
    pal.table <- palindrome[,c(1,3,4)]
    mpal.table <- m.palindrome[,c(1,3,4)]
    names(pal.table) <- NULL
    names(mpal.table) <- NULL
    write.table(pal.table,"all.unmasked.pal.tab",sep="\t",
                row.names = FALSE,quote = FALSE)
    write.table(mpal.table,"all.masked.pal.tab",sep="\t",
                row.names = FALSE,quote = FALSE)
