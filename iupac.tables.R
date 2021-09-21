#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script

source("/home/sq16564/Analysis-of-pig-X/Rfunctions.R") #Call functions pre-defined in script 

palindrome.tab <- read.table(args[1], header = FALSE, sep = "") #Read table 1 with palindrome start and stop positions for the first half of the unmasked DNA sequence
palindrome.tab2 <- read.table(args[2], header = FALSE, sep = "") #Read table 2 with palindrome start and stop positions for the second half of the unmasked DNA sequence
palindrome.tab <- palindrome.tab %>% dplyr::rename(start = V1,stop = V2,start2 = V3,stop2 = V4) #Rename headers of table to correspond to start and stop positions 1 and 2 (applied for later use of table in scripts)
palindrome.tab2 <- palindrome.tab2 %>% dplyr::rename(start = V1,stop = V2,start2 = V3,stop2 = V4) #Rename headers of table to correspond to start and stop positions 1 and 2 (applied for later use of table in scripts)
palindrome.table <- rbind(palindrome.tab,palindrome.tab2) #Combine the two tables of palindrome sequences for the two halfs of the DNA sequence 
write.table(palindrome.tab, file = paste0(args[1]), quote = FALSE, sep ="\t",row.names = FALSE) #Write renamed table 
