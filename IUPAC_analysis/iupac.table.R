#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script

source("/home/sq16564/Analysis-of-pig-X/Rfunctions.R") #Call functions pre-prepared in Rfunctions script 

palindrome.tab <- read.table(args[1], header = FALSE, sep = "") #Load palindrome position table into R 
palindrome.tab <- palindrome.tab %>% dplyr::rename(start = V1,stop = V2,start2 = V3,stop2 = V4) #Rename headers of table to correspond to start and stop positions 1 and 2 (applied for later use of table in scripts)
write.table(palindrome.tab, file = paste0(args[1]), quote = FALSE, sep ="\t",row.names = FALSE) #Write renamed table 
