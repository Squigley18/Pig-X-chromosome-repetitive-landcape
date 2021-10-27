#Taking the subject start and stop positions from the BLAST results tables they can be input into these functions to create overlap line graphs... 
#individual functions were created for each subject accession region to allow for annotations to be added corresponding to the NCBI annotations


overlaps.to.L1 <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
L1.acc <- make.windows(1,7878)
overlap.hits <- countOverlaps(L1.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:7878)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
transposase <- as.data.frame(cbind(start=2324,stop=3214,yaxis=-5))
reverse.transcriptase <- as.data.frame(cbind(start=3740,stop=7558,yaxis=-5))
repeat.1 <- as.data.frame(cbind(start=1230,stop=2080,yaxis=-5))
repeat.2 <- as.data.frame(cbind(start=400,stop=680,yaxis=-5))
UTR5 <- as.data.frame(cbind(start=15,stop=2330,yaxis=-11))
UTR3 <- as.data.frame(cbind(start=7595,stop=7865,yaxis=-11))
LINE <- as.data.frame(cbind(start=15,stop=7865,yaxis=-11))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=LINE,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkblue", col=NA, inherit.aes=F)+
	geom_rect(data=transposase,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=reverse.transcriptase,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=UTR5,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="grey39", col=NA, inherit.aes=F)+
	geom_rect(data=UTR3,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="grey39", col=NA, inherit.aes=F)+
	geom_rect(data=repeat.1,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="mediumblue", col=NA, inherit.aes=F)+
	geom_rect(data=repeat.2,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="mediumblue", col=NA, inherit.aes=F)+
	scale_x_continuous(name="L1 accession region (bp)",label=c(0,2500,5000,7500), breaks=c(0,2500,5000,7500))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(3940,7730,1172,540,1655,5649,2769), y = c(-7,-7,-7,-1,-1,-1,-1),
             label = c("L1", "3'UTR", "5'UTR","repeat","repeat","reverse transcriptase", "transposase") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.L1("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 



overlaps.to.WIF1 <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,401028)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:401028)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
WIF.gene <- as.data.frame(cbind(start=1,stop=61511,yaxis=-10))
LEM.gene <- as.data.frame(cbind(start=113638,stop=195928,yaxis=-10))
MSRB3.gene <- as.data.frame(cbind(start=218130,stop=386195,yaxis=-10))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=WIF.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=LEM.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=MSRB3.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	scale_x_continuous(name="WIF1 accession region (kbp)",label=c(0,100,200,300,400), breaks=c(0,100000,200000,300000,400000))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(30755,154783,302162), y = c(-6,-6,-6),
             label = c("WIF1", "LEMD3", "MSRB3") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.WIF1("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.to.myostatin <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,167404)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:167404)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
myo.gene <- as.data.frame(cbind(start=81992,stop=86906,yaxis=-5))
myo.exon <- as.data.frame(cbind(start=81992,stop=82364,yaxis=-5))
myo.exon2 <- as.data.frame(cbind(start=84174,stop=84547,yaxis=-5))
myo.exon3 <- as.data.frame(cbind(start=86526,stop=86906,yaxis=-5))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=myo.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=myo.exon,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	geom_rect(data=myo.exon2,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	geom_rect(data=myo.exon3,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Myostatin accession region (bp)",label=c(0,41851,83702,125553,167404), breaks=c(0,41851,83702,125553,167404))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(84449), y = c(-1),
             label = c("myostatin") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.myostatin("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.to.pig.histone <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,2646)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:2646)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
histone.loci <- as.data.frame(cbind(start=1,stop=2646,yaxis=-10))
histone.gene <- as.data.frame(cbind(start=1536,stop=1877,yaxis=-5))
histone.exon <- as.data.frame(cbind(start=1556,stop=1646,yaxis=-5))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=histone.loci,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	geom_rect(data=histone.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=histone.exon,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Pig histone H2A (bp)",label=c(0,662,1323,1984,2646), breaks=c(0,662,1323,1984,2646))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(1323,1706), y = c(-6,-1),
             label = c("LOC100517800","Histone H2A type 2/3 like") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.pig.histone("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.to.lagen.histone <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,552)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:552)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
lagen.loci <- as.data.frame(cbind(start=1,stop=552,yaxis=-11))
lagen.gene <- as.data.frame(cbind(start=1,stop=552,yaxis=-6))
lagen.exon <- as.data.frame(cbind(start=16,stop=115,yaxis=-6))


  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=lagen.loci,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	geom_rect(data=lagen.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=lagen.exon,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Lagenorhynchus obliquidens histone H2A (bp)",label=c(0,138,276,414,552), breaks=c(0,138,276,414,552))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(277,277), y = c(-7,-2),
             label = c("LOC113616392","Histone H2A type 2/3 like") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.lagen.histone("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.pig.full.loci <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,1746)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:1746)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
pig.loci <- as.data.frame(cbind(start=1,stop=1746,yaxis=-6))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=pig.loci,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Pig uncharacterised loci LOC110257707 (bp)",label=c(0,437,873,1309,1746), breaks=c(0,437,873,1309,1746))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(873), y = c(-2),
             label = c("LOC110257707") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.pig.full.loci("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.camel.full.loci <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,917)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:917)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
camel.loci <- as.data.frame(cbind(start=1,stop=917,yaxis=-6))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=camel.loci,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Camel uncharacterised loci LOC116152181 (bp)",label=c(0,230,459,687,917), breaks=c(0,230,459,687,917))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(459), y = c(-2),
             label = c("LOC116152181") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.camel.full.loci("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.DCT <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,127098)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:127098)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
DCT.gene <- as.data.frame(cbind(start=13585,stop=60979,yaxis=-6))
GPC6.gene <- as.data.frame(cbind(start=91652,stop=111903,yaxis=-6))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=DCT.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=GPC6.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	scale_x_continuous(name="DCT accession region (bp)",label=c(0,31750,63500,95250,127000), breaks=c(0,31750,63500,95250,127000))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(37282,101777), y = c(-2,-2),
             label = c("DCT","GPC6") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.DCT("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.CBX6 <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,163158)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:163158)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
CBX6.gene <- as.data.frame(cbind(start=61860,stop=67422,yaxis=-6))
APOBEC3.gene <- as.data.frame(cbind(start=132076,stop=149427,yaxis=-6))
APOBEC3.exons <- as.data.frame(cbind(start=c(132067,134520,138167,138855,139461,139679,142439,147172,148256,148772,149211),
stop=c(132379,134670,138446,138969,139574,140037,142586,147328,148523,148887,149427),yaxis=-6))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=CBX6.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=APOBEC3.gene,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=APOBEC3.exons,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	scale_x_continuous(name="CBX6 accession region (bp)",label=c(0,40750,81500,122250,163000), breaks=c(0,40750,81500,122250,163000))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(64641,140751), y = c(-2,-2),
             label = c("CBX6","APOBEC3") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.CBX6("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.pervA <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,4561)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:4561)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
LTR1 <- as.data.frame(cbind(start=1248,stop=2021,yaxis=-6))
LTR2 <- as.data.frame(cbind(start=2122,stop=2754,yaxis=-6))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=LTR1,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkblue", col=NA, inherit.aes=F)+
	geom_rect(data=LTR2,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkblue", col=NA, inherit.aes=F)+
	scale_x_continuous(name="PERV-A accession region (bp)",label=c(0,1140,2280,3420,4561), breaks=c(0,1140,2280,3420,4561))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(1634,2438), y = c(-2,-2),
             label = c("LTR","LTR") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.PERVA("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.pervC <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,8955)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:8955)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
LTR1 <- as.data.frame(cbind(start=282,stop=893,yaxis=-6))
LTR2 <- as.data.frame(cbind(start=8344,stop=8955,yaxis=-6))
gag <- as.data.frame(cbind(start=1343,stop=2917,yaxis=-11))
pro.pol <- as.data.frame(cbind(start=3065,stop=6502,yaxis=-6))
env <- as.data.frame(cbind(start=6378,stop=8294,yaxis=-11))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=LTR1,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkblue", col=NA, inherit.aes=F)+
	geom_rect(data=LTR2,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkblue", col=NA, inherit.aes=F)+
	geom_rect(data=gag,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=pro.pol,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=env,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	scale_x_continuous(name="PERV-C accession region (bp)",label=c(0,2225,4450,6675,8900), breaks=c(0,2225,4450,6675,8900))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(587,8649,4788,2157,7336), y = c(-2,-2,-2,-7,-7),
             label = c("LTR","LTR","pro-pol","gag","env") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.pervC("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.pig.staggered.loci <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,3035)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:3035)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
pig.loci <- as.data.frame(cbind(start=1,stop=3035,yaxis=-6))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=pig.loci,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Pig uncharacterised loci LOC110257827 (bp)",label=c(0,750,1500,2250,3035), breaks=c(0,750,1500,2250,3035))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(1518), y = c(-1),
             label = c("LOC110257827") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.pig.staggered.loci("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.cow.factor <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,4916)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:4916)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
cow.F8A1 <- as.data.frame(cbind(start=1,stop=4916,yaxis=-7))
cow.intron22.protein <- as.data.frame(cbind(start=56,stop=1138,yaxis=-7))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=cow.F8A1,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=cow.intron22.protein,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="yellow2", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Cattle coagulation factor 8 associated 1 (bp)",label=c(0,1225,2490,3675,4916), breaks=c(0,1225,2490,3675,4916))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(2458,597), y = c(-2,-2),
             label = c("F8A1","inron 22 protein") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.cow.factor("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 


overlaps.cat.factor <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
wif.acc <- make.windows(1,1274)
overlap.hits <- countOverlaps(wif.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:1274)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
cat.loci <- as.data.frame(cbind(start=1,stop=1274,yaxis=-20))
cat.intron22.protein <- as.data.frame(cbind(start=46,stop=1173,yaxis=-12))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=cat.loci,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	geom_rect(data=cat.intron22.protein,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="yellow2", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Cat factor 8 intron 22 protein (bp)",label=c(0,319,637,956,1274), breaks=c(0,319,637,956,1274))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(637,610), y = c(-15,-6),
             label = c("LOC101095239","intron 22 protein") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.cat.factor("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.to.HSFX <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
HSFX.acc <- make.windows(1,1216)
overlap.hits <- countOverlaps(HSFX.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:1216)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
locus <- as.data.frame(cbind(start=1,stop=1216,yaxis=-10))
gene.HSFX <- as.data.frame(cbind(start=218,stop=1138,yaxis=-5))
exon <- as.data.frame(cbind(start=417,stop=486,yaxis=-5))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=locus,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	geom_rect(data=gene.HSFX,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=exon,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Ferret HSFX accession region (bp)",label=c(0,304,608,912,1216), breaks=c(0,304,608,912,1216))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(678,609), y = c(-1,-7),
             label = c("HSFX","LOC106003320") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.HSFX("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 

overlaps.to.camel.HSFX <- function(tab){
subject <-  read.table(tab, header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) 

subject <- subject %>% dplyr::rename(start = V1,stop = V2) 
chromosome <- "X"
strand <- "*" 
hits <- cbind(subject,chromosome,strand)
hitsGF <- GRfromDF(hits)
HSFX.acc <- make.windows(1,4698)
overlap.hits <- countOverlaps(HSFX.acc,hitsGF)
overlaps <- as.data.frame(overlap.hits)
xaxis <- seq(1:4698)
overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
locus <- as.data.frame(cbind(start=1,stop=4698,yaxis=-10))
gene.HSFX <- as.data.frame(cbind(start=1655,stop=2899,yaxis=-5))
exon <- as.data.frame(cbind(start=1747,stop=1861,yaxis=-5))

  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
	geom_rect(data=locus,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkgreen", col=NA, inherit.aes=F)+
	geom_rect(data=gene.HSFX,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
	geom_rect(data=exon,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="black", col=NA, inherit.aes=F)+
	scale_x_continuous(name="Camel HSFX accession region (bp)",label=c(0,1174,2348,3522,4698), breaks=c(0,1174,2348,3522,4698))+
    theme_classic()+
    ylab("Number of overlapping LASTZ hits")+
	annotate("text", x = c(2277,2349), y = c(-1,-7),
             label = c("HSFX","LOC106003320") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

plot <- overlaps.to.camel.HSFX("data") #input subject start and stop positions 
ggsave("data.svg",plot,unit="mm",height=85,width=120) # export as SVG graph 
