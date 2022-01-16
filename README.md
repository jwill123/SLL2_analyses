# SLL2_analyses
Code for analyzing the data of the 2nd edition of Saca La Lengua 

The primary scripts for the analyses are:

- **src/SLL2_dada2.R** for read filtering and taxonomy assignment
- **src/SLL_2_phyloseq.R** for statistical analyses and figure generation
- **src/SLL_2_phyloseq.functions.R** for the code for the many of the functions used in src/SLL_2_phyloseq.R

The raw fastq files can be found in the Sequence Read Archive (SRA) with the accession number **PRJNA667146** and can be found here: http://www.ncbi.nlm.nih.gov/bioproject/667146.

The additional data necessary for the analyses are:

- **data/SLL2.metadata.xlsx** - a table of the metadata used for the analyses
- **data/SLL2.MALDI_results.xlsx** - the results of MALDI-TOF analyses for the identification of yeast species in our samples
