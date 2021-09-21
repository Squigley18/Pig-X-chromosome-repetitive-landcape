# Pig-X-chromosome-repetitive-landcape
Am I repeating myself? Determining the repetitive landscape of the pig X chromosome. 

# Introduction
This git repository contains all the scripts used in the analysis of the pig X chromosome. This study aimed to determine the repetitive content of the pig X chromosome involving finding inverted repeats (IUPACpal analysis) and finding regions of high similarity in a self-alignment (LASTZ alignment and analysis). 
The pig DNA sequence investigated was assessed both with repeatmasking applied and the unmasked sequence. The LASTZ alignments were also filtered to 95% and 99% identity between alignments. 
These alignments were then assessed further using R scripts to determine regions of the chromsome with high numbers of alignment hits. This then led to extracting all the known genes in the X chromosome and determining which of these genes had overlapping alignment hits. These alignment hits overlapping genes of the X chromosome were plot to determine their distribution. Finally the start and stop positions of these alignment hits were exported to a table which was used to extract FASTA sequences from the pig X chromosome to run through the ncbi BLAST software.
Running the DNA sequences through BLAST determined genes which showed similarity to these sequences of DNA. This then provided a table of alignment hits with the DNA sequence within the IUPAC or LASTZ hits and the regions of similarity within the subject sequence. This results table allows for the exact regions of the chromosome where this similarity is found to be extracted in FASTA format and the DNA sequence can be manually aligned to the homologous genes to determine the extent of the similarity and what genomic elements are found within those regions. 
Combining these results the regions of self homology in the X chromosome can be investigated to understand their nature, their distribution which allows for determining the nature of the replication, and their significance particularly where the regions of homology are found in introns or exons. Here, any amplicomic genes, repetitive elements or duplicated genes can be identified. 

# Folders in the repository 

| FILE / FOLDER | PURPOSE |
| :--- | :--- |
| `functions` | Folder of functions created to be used within scripts 
| `IUPAC_analysis` | Folder of scripts used for the analysis of inverted repeats in pig X chromosome 
| `LASTZ_alignment` | Folder of scripts used for self-alignment of pig X chromosome
| `LASTZ_analysis` | Folder of scripts used in computational analysis of LASTZ data
| `BLAST_analysis` | Folder of scripts used in the ncbi BLAST analysis of pig DNA sequences

# Functions 
The following scripts contain all the functions created for the analysis of the IUPAC data and LASTZ data in the self-alignment of the pig X chromosome

| SCRIPT | OVERVIEW |
| :--- | :--- | 
| `Rfunctions.R` | Functions created for the R scripts in the analysis of the LASTZ data 
| `Rfunctions.iupac.R` | Functions created for the R scripts in the analysis of the IUPAC data similar to the LASTZ functions 
| `shell-script-functions.sh` | Functions created for the shell scripts in the analysis of the LASTZ and IUPAC data

Some of the functions include: 

```
overlap.hits.to.chromosome()
```
Overlap.hits.to.chromosome function is used to create the dotplot showing the number of alignment hits sound in regions of the X chromosome 

```
subset.overlap.segment.plot()
```
Subset.overlap.segment.plot function was created to plot the genes of the X chromosome as a segment in red with the alignment hits overlapping the gene in black. The gene then is annotated to highlight regions of inrons, exons and untranslated regions. 

```
Ideogram.plot()
``` 
Ideogram.plot function creates an outline of the X chromosome ideogram, and is then overlapped with vertical segments showing the distribution of the alignment hits. 

# IUPAC analysis

| SCRIPTS | OVERVIEW |
| :--- | :--- | 
| `split.dna.files.sh` | Scripts to split the pig DNA sequence in two before the IUPAC analysis to reduce computational demand
| `iupacpal.sh` |  script for IUPAC analysis producing a text based format output file 
| `iupac.table.sh` | Script to convert text based format to table which can be analysed in R scripts
| `iupac.table.masked.sh` | Script to convert text based format to table which can be analysed in R scripts
| `iupac.table.R` | R script to add a header to the table which will be consistent with analysis script
| `iupac.tables.R` | R script to add a header to the table which will be consistent with analysis script and combine the unmasked DNA tables into one
| `overlap.hits.entire.chr.iupac.R` | R script to arrange the IUPAC palindromes on a graph to show the number of palindromes found in each region of the chromosome
| `entire.chromosome.iupac.R` | R script to make overlapping segment graphs between genes and palindrome hits, the distribution of the hits along an ideogram and exporting a table of hit coordinates 
| `iupac.blast.filtering.sh` |  Script to run the palindrome hits through ncbi BLAST to find similar genes and then extract the specific regions of similarity from the X chromosome FASTA sequence for further analysis
| `perl.extract.tsv.files.pl` | Perl script to extract DNA sequence from pig X chromosome using Ensembl API using given coordinates
| `wrapper.iupac.sh` | Combination of the above scripts within one wrapper script to be used in one 'smooth' run

The above scripts were developed to optimise the analysis of the inverted repeat content of the pig X chromosome. The DNA sequence files were analysed with repeatmasking applied and unmasked also. These files are around 130MB in size and lead to high computational demand, therefore they were split in half and run as seperate files. The IUPAC output is returned as a text based format therefore many of the scripts involve reformatting this to allow the data to be analysed further. The detected palindromes are then plot to determine the number of hits in the X chromosome and find any genes which these hits overlap. The complimentary palindromes were also plot to determine their distribution throughout the chromosome. Finally the DNA sequence within these palindromes is extracted to be run through the ncbi BLAST software to determine their nature. 

# LASTZ alignment 

| SCRIPTS | OVERVIEW |
| :--- | :--- | 
| `LASTZ-unmasked.alignment.sh` | Script for self alignment of pig X chromosome using the unmasked DNA sequence producing an alignment table
| `LASTZ-masked.alignment.sh` | Script for self alignment of pig X chromosome using the masked DNA sequence producing an alignment table 
| `LASTZ-unmasked.dotplot.sh` | Script for self alignment of pig X chromosome using the unmasked DNA sequence producing data to create a naieve dot plot
| `LASTZ-masked.dotplot.sh` | Script for self alignment of pig X chromosome using the masked DNA sequence producing data to create a naieve dot plot
| `LASTZ-dotplot.R` | R script to produce the naieve dot plot png graph 
| `naieve.r.dotplot.sh` | wrapper script to produce the naieve dot plots for the masked and unmasked DNA sequences at the different percentage identities 

The above scripts were written to utilise the LASTZ function to perform the LASTZ self alignment of the pig DNA sequence. The alignment performed set the percentage identity of the LASTZ hits to be above 95% and 99% identity. The alignments of the masked and unmasked DNA sequences required different scripts due to the unmasked DNA sequence requiring the extra hspthresh flag to combat the high computational demand of the alignment. The initial alignment files provided a table of coordinates of the LASTZ alignment hits and the dotplot scripts produce a table of hits which can be used to create a naieve dot matrix plot. The alignment table can then be used in further analysis processes. 


