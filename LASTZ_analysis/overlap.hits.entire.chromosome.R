#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 

source("/home/sq16564/Analysis-of-pig-X/Rfunctions.R") #call script with R functions created 
  
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
  
#Count the overlaps between the LASTZ genomic ranges objects and the chromosome window of 3000bp in length (chosen due to average LASTZ hit length)
overlap.ID99 <- countOverlaps(window.3000,ID99gr)
overlap.ID95 <- countOverlaps(window.3000,ID95gr)
overlap.IDrm99 <- countOverlaps(window.3000,IDrm99gr)
overlap.IDrm95 <- countOverlaps(window.3000,IDrm95gr)

# The X chromosome GenomicRanges object has 42000 rows which represents the windows of 3000bp along the X chromosome 

xaxis <- seq(1:42000)

# Convert the countOverlaps objects into data frames which can then be used to create the dotplot through ggplot where the X axis shows the location on the X chromosome and the Y axis the count for the number of overlapping LASTZ hits in that region 
# Any regions with 0 overlaps are converted to NA in the dataframe to allow for removing the NA values from the ggplot

overlap.ID99.df <- as.data.frame(cbind(xaxis,overlap.ID99))
overlap.ID99.df[overlap.ID99.df==0] <- NA
overlap.ID99.df <- overlap.ID99.df %>% dplyr::rename(yaxis = overlap.ID99)

overlap.ID95.df <- as.data.frame(cbind(xaxis,overlap.ID95))
overlap.ID95.df[overlap.ID95.df==0] <- NA
overlap.ID95.df <- overlap.ID95.df %>% dplyr::rename(yaxis = overlap.ID95)

overlap.IDrm99.df <- as.data.frame(cbind(xaxis,overlap.IDrm99))
overlap.IDrm99.df[overlap.IDrm99.df==0] <- NA
overlap.IDrm99.df <- overlap.IDrm99.df %>% dplyr::rename(yaxis = overlap.IDrm99)


overlap.IDrm95.df <- as.data.frame(cbind(xaxis,overlap.IDrm95))
overlap.IDrm95.df[overlap.IDrm95.df==0] <- NA
overlap.IDrm95.df <- overlap.IDrm95.df %>% dplyr::rename(yaxis = overlap.IDrm95)

# The following are the objects for the dot plots created with the overlap.hits.to.chromosome function which can then be exported as pngs using ggsave 

ID99.count.plot <- overlap.hits.to.chromosome(overlap.ID99.df,"99% identity hits distribution")
ID95.count.plot <- overlap.hits.to.chromosome(overlap.ID95.df,"95% identity hits distribution")
IDrm99.count.plot <- overlap.hits.to.chromosome(overlap.IDrm99.df,"repeat-masked 99% identity hits distribution")
IDrm95.count.plot <- overlap.hits.to.chromosome(overlap.IDrm95.df,"repeat-masked 95% identity hits distribution")

ggsave(paste0(args[1],".overlap.svg"),ID99.count.plot,dpi=300,unit="mm",height=85,width=120)
ggsave(paste0(args[2],".overlap.svg"),ID95.count.plot,dpi=300,unit="mm",height=85,width=120)
ggsave(paste0(args[3],".overlap.svg"),IDrm99.count.plot,dpi=300,unit="mm",height=85,width=120)
ggsave(paste0(args[4],".overlap.svg"),IDrm95.count.plot,dpi=300,unit="mm",height=85,width=120)

# manipulate the LASTZ data to calculate the overlapping hits in regions of the chromosome
I99.hits <- manipulate.data.for.overlap.ideogram.plot(ID99)
Rm99.hits <- manipulate.data.for.overlap.ideogram.plot(IDrm99)
I95.hits <- manipulate.data.for.overlap.ideogram.plot(ID95)
Rm95.hits <- manipulate.data.for.overlap.ideogram.plot(IDrm95)
#Combining the hit objects into one dataframe to be passed to the function to create the ideogram plot 
all.ideogram.hits <- as.data.frame(cbind(I99.hits$xaxis,I99.hits$ideogram,Rm99.hits$ideogram,I95.hits$ideogram,Rm95.hits$ideogram))
all.ideogram.hits <- all.ideogram.hits %>% dplyr::rename(Xaxis = V1, ID99 = V2, Rm99 = V3, ID95 = V4, Rm95 = V5)
#Create and save the ideogram plot 
all.plot <- Ideogram.overlap.plot(all.ideogram.hits,"all palindrome hits distribution")
ggsave("all.LASTZ.hits.ideogram.svg",all.plot,unit="mm",height=85,width=120)
