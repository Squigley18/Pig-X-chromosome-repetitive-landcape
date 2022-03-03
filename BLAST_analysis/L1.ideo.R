#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script 
  
source("/home/sq16564/Analysis-of-pig-X/Rfunctions.R")
  
 # edit manipulate data function to work with L1 data  
  
manipulate.data.for.overlap.ideogram.plot <- function(mydf){
  hit.gr <- GRfromDF(mydf)
  window.3000 <- make.windows(3000,126e6)
  hits <- countOverlaps(window.3000,hit.gr)
  xaxis <- seq(1:42000)
  counts <- as.data.frame(cbind(xaxis,hits))
  counts$ideogram <- counts$hits
  counts$ideogram[counts$hits >= 1] <- counts$xaxis[counts$hits >= 1]
  counts$ideogram[counts$ideogram==0] <- NA
  counts
}


  
  chromosome.ideogram <- data.frame(x = 1:42000, y = c(rep.int(1,20500),
                                                       rep.int(NA,800),rep.int(1,20700)),
                                    z = c(rep.int(2,20500),
                                          rep.int(NA,800),rep.int(2,20700)))
  
  
Rm99.hits <- read.table(args[1])
Rm95.hits <- read.table(args[2])

colnames(Rm99.hits) <- c("chromosome","start","stop")
colnames(Rm95.hits) <- c("chromosome","start","stop")

ideo.99.L1 <- manipulate.data.for.overlap.ideogram.plot(Rm99.hits)
ideo.95.L1 <- manipulate.data.for.overlap.ideogram.plot(Rm95.hits)

Ideogram <- ggplot(chromosome.ideogram, aes(x,y),na.rm = TRUE)+ 
  geom_path(size = 4, lineend = "round",colour="gray87")+
  geom_path(aes(x,z),size = 4, lineend = "round",colour="gray87")+
  geom_segment(data = ideo.95.L1,aes(x=ideogram,xend=ideogram+25,
                                     y=1,yend=1),size=3.5,colour="#1026EB")+
  geom_segment(data = ideo.99.L1,aes(x=ideogram,xend=ideogram+25,
                                     y=2,yend=2),size=3.5,colour="#0A0A0A")+
  scale_x_continuous(name="Location on chromosome (Mbp)", label=c(0,30,60,90,120), breaks=c(0,10000,20000,30000,40000))+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("L1_ideogram.svg",Ideogram,unit="mm",height=85,width=120)
