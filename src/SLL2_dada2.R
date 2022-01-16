


library(dada2)
library(ggplot2)


# depending on where the script is located...
if (dir.exists("/users/tg/jwillis/SLL")) {
  tg_dir   <- "/users/tg"
  seq_dir_2018  <- "/users/tg/sequencing_data/ESaus/"
  seq_dir_2019  <- "/users/tg/sequencing_data/ESaus/"
  home_dir <- "/users/tg/jwillis/SLL"
} else if (dir.exists("/gpfs/projects/bsc40/current/jwillis/SLL")) {
  tg_dir   <- "/gpfs/projects/bsc40"
  seq_dir_2018  <- "/gpfs/projects/bsc40/sequencing_data/repository/Microbioma/2018"
  seq_dir_2019  <- "/gpfs/projects/bsc40/sequencing_data/repository/Microbioma/2019"
  home_dir <- "/gpfs/projects/bsc40/current/jwillis/SLL"
} else if (dir.exists("/home/jwillis/gpfs/projects/bsc40/current/jwillis/SLL")) {
  tg_dir   <- "/home/jwillis/gpfs/projects/bsc40"
  seq_dir_2018  <- "/home/jwillis/gpfs/projects/bsc40/sequencing_data/repository/Microbioma/2018"
  seq_dir_2019  <- "/home/jwillis/gpfs/projects/bsc40/sequencing_data/repository/Microbioma/2019"
  home_dir <- "/home/jwillis/gpfs/projects/bsc40/current/jwillis/SLL"
} else if (dir.exists(sprintf("%s/SLL", getwd()))) {
  tg_dir   <- sprintf("%s/../", getwd())
  seq_dir  <- ""
  home_dir <- sprintf("%s/SLL", getwd())
} else if (dir.exists("~/Downloads/SLL")) {
  tg_dir   <- ""
  seq_dir  <- ""
  home_dir <- "~/Downloads/SLL"
} else {
  tg_dir   <- ""
  seq_dir  <- ""
  home_dir <- getwd()
}


### SEE TUTORIAL HERE:   https://benjjneb.github.io/dada2/tutorial.html


# ****************************************************** #
path <- sprintf("%s/Part_2/reads", home_dir)
length(list.files(path))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: datax_xxx_SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



# ****************************************************** #
#### EXAMINE QUALITY PROFILES OF FORWARD AND REVERSE READS #### 
# ****************************************************** #
plotQualityProfile(fnFs[1:4]) + 
  scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5))
# From this it seems should truncate forward reads at least at position 250 (trim last 50 nucs), and trim first 10

plotQualityProfile(fnRs[1:4]) +
  scale_x_continuous(breaks = seq(0, 300, 25)) +
  scale_y_continuous(breaks = seq(0, 40, 5))
# reverse have worse quality at ends, truncate at position 225, and trim first 10




# ****************************************************** #
#### PERFORM FILTERING AND TRIMMING #### 
# ****************************************************** #

# This removed almost half the reads per sample on average
# So here I output filtered reads into a different folder, using a more relaxed maxEE filter for the reverse reads
filt_path <- file.path(sprintf("%s/Part_2/reads", home_dir), "filtered_reads") # Place filtered files in filtered_reads/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(250,225), maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                     trimLeft=c(10,10), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

mean(out[,2] / out[,1]) # 0.6423552 -- not bad


# *********************************** #


filt_path.2 <- file.path(sprintf("%s/Part_2/reads", home_dir), "filtered_reads_2") # Place filtered files in filtered_reads/ subdirectory
filtFs.2 <- file.path(filt_path.2, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs.2 <- file.path(filt_path.2, paste0(sample.names, "_R_filt.fastq.gz"))

out.2 <- filterAndTrim(fnFs, filtFs.2, fnRs, filtRs.2,
                     truncLen=c(275,230), maxN=0, maxEE=c(10,10), truncQ=2, rm.phix=TRUE,
                     trimLeft=c(10,10), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

mean(out.2[,2] / out.2[,1]) # 0.7253076 -- not bad




# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))



# ****************************************************** #
#### LEARN THE ERROR RATES #### 
# ****************************************************** #
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

errF.2 <- learnErrors(filtFs.2, multithread=TRUE)
errR.2 <- learnErrors(filtRs.2, multithread=TRUE)

# when this is run with the filtFs, stops after sample 149 because the 'nbases' parameter is 1e8 reads
#   by default, stops after this many are read into memory. Obviously there is a lot of data remaining,
#   so if the values look weird after running this, maybe should increase the 'nbases' value


# show error rates for each possible transition (eg. A->C, A->G, …) 
plotErrors(errF, nominalQ=TRUE) # looks pretty good
plotErrors(errR, nominalQ=TRUE) # looks all good

plotErrors(errF.2, nominalQ=TRUE) # looks pretty good
plotErrors(errR.2, nominalQ=TRUE) # looks all good


# Points are the observed error rates for each consensus quality score. 
# The black line shows the estimated error rates after convergence. 
# The red line shows the error rates expected under the nominal definition of the Q-value.

# **** If using this workflow on your own data ****:
# Parameter learning is computationally intensive, so by default the learnErrors function 
#    uses only a subset of the data (the first 1M reads). If the plotted error model does 
#    not look like a good fit, try increasing the nreads parameter to see if the fit improves.



# from this post about not converging after 10 rounds or error estimation https://github.com/benjjneb/dada2/issues/77
dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)

dada2:::checkConvergence(errF.2)
dada2:::checkConvergence(errR.2)





# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))





# ****************************************************** #
#### DEREPLICATION, SAMPLE INFERENCE, MERGE PAIRED READS #### 
# ****************************************************** #
# Combines all identical sequences into unique sequences with abundances
# Reduces computation time by eliminating redundant comparisons
# *** UNIQUE TO DADA2: retains summary of quality info for each unique sequence
#     This is an average of positional qualities for duplicated reads to inform the error model

# Since there are so many samples, better to process them in streaming fashion (for loop)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  # cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[ grep(sam, filtFs) ], verbose = T)
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  
  derepR <- derepFastq(filtRs[ grep(sam, filtRs) ], verbose = T)
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
# rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, sprintf("%s/Part_2/DADA2/seqtab.rds", home_dir))


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# only keep 
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(298,310)]


# *********************************** #

mergers.2 <- vector("list", length(sample.names))
names(mergers.2) <- sample.names

# for(sam in sample.names[1:20]) {
# for(sam in sample.names[21:60]) {
# for(sam in sample.names[61:100]) {
# for(sam in sample.names[101:150]) {
# for(sam in sample.names[151:200]) {
# for(sam in sample.names[201:280]) {
# for(sam in sample.names[281:340]) {
# for(sam in sample.names[341:400]) {
# for(sam in sample.names[401:450]) {
# for(sam in sample.names[451:520]) {
# for(sam in sample.names[521:600]) {
# for(sam in sample.names[601:length(sample.names)]) {
for(sam in sample.names) {
  # cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs.2[ grep(sam, filtFs.2) ], verbose = T)
  ddF <- dada(derepF, err=errF.2, multithread=TRUE)
  
  derepR <- derepFastq(filtRs.2[ grep(sam, filtRs.2) ], verbose = T)
  ddR <- dada(derepR, err=errR.2, multithread=TRUE)
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers.2[[sam]] <- merger
}
# rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab.2 <- makeSequenceTable(mergers.2)
saveRDS(seqtab.2, sprintf("%s/Part_2/DADA2/seqtab.2.rds", home_dir))


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.2)))




# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))







# ****************************************************** #
#### MERGE SEQUENCE RUNS (if necessary), REMOVE CHIMERAS, ASSIGN TAXONOMY #### 
# ****************************************************** #
# This takes the seqtab (or the merged seqtabs from the different sequencing runs) and removes chimeras
# Then assign taxonomy as usual

# Merge multiple runs (if necessary)
# st1 <- readRDS("/users/tg/jwillis/SLL/Part_1/DADA2/seqtab1.rds")
# st2 <- readRDS("/users/tg/jwillis/SLL/Part_1/DADA2/seqtab2.rds")
# st3 <- readRDS("/users/tg/jwillis/SLL/Part_1/DADA2/seqtab3.rds")
# st.all <- mergeSequenceTables(st1, st2, st3)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
dim(seqtab)
dim(seqtab.nochim)
ncol(seqtab.nochim)/ncol(seqtab)
1 - sum(seqtab.nochim)/sum(seqtab)

# Shows that 3743 / 5490 seqs (68.17%) were removed as bimeras, but they only make up ~2.3% of total seq reads


# Assign taxonomy
path.tut <- "/users/tg/jwillis/SLL/DADA2_files"
tax <- assignTaxonomy(seqtab.nochim, sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)

# Add species
tax_spec <- addSpecies(tax, sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))

# remove "." from any taxa names because it will break part of the pipeline later that relies on splitting names by "."
tax[ grepl('\\.',tax) ] <- gsub('\\.','_',tax[grepl('\\.',tax)])
tax_spec[ grepl('\\.',tax_spec) ] <- gsub('\\.','_',tax_spec[grepl('\\.',tax_spec)])

# Write to disk
saveRDS(seqtab.nochim, "/users/tg/jwillis/SLL/Part_2/DADA2/seqtab_final.rds")
saveRDS(tax, "/users/tg/jwillis/SLL/Part_2/DADA2/tax_final.rds")
saveRDS(tax_spec, "/users/tg/jwillis/SLL/Part_2/DADA2/tax_species_final.rds")


# inspect taxonomic assignments
tax.print <- tax # Removing sequence rownames for display only
rownames(tax.print) <- NULL
head(tax.print)

tax_spec.print <- tax_spec # Removing sequence rownames for display only
rownames(tax_spec.print) <- NULL
head(tax_spec.print)

summary(is.na(tax_spec.print))


# *********************************** #


# Remove chimeras
seqtab.nochim.2 <- removeBimeraDenovo(seqtab.2, method="consensus", multithread=TRUE)
dim(seqtab.2)
dim(seqtab.nochim.2)
ncol(seqtab.nochim.2)/ncol(seqtab.2)
1 - sum(seqtab.nochim.2)/sum(seqtab.2)

# Shows that 10526 / 32983 seqs (31.9%) were removed as bimeras, but they only make up ~16% of total seq reads


# Assign taxonomy
path.tut <- "/users/tg/jwillis/SLL/DADA2_files"
tax.2 <- assignTaxonomy(seqtab.nochim.2, sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)

# Add species
tax_spec.2 <- addSpecies(tax.2, sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))
# change '.' to '_' in any names to avoid problems in later parts
tax_spec.2[ grepl('\\.',tax_spec.2) ] <- gsub('\\.','_',tax_spec.2[grepl('\\.',tax_spec.2)])

# Write to disk
saveRDS(seqtab.nochim.2, "/users/tg/jwillis/SLL/Part_2/DADA2/seqtab_final.2.rds")
saveRDS(tax.2, "/users/tg/jwillis/SLL/Part_2/DADA2/tax_final.2.rds")
saveRDS(tax_spec.2, "/users/tg/jwillis/SLL/Part_2/DADA2/tax_species_final.2.rds")


# inspect taxonomic assignments
tax.print.2 <- tax.2 # Removing sequence rownames for display only
rownames(tax.print.2) <- NULL
head(tax.print.2)

tax_spec.print.2 <- tax_spec.2 # Removing sequence rownames for display only
rownames(tax_spec.print.2) <- NULL
head(tax_spec.print.2)

summary(is.na(tax_spec.print.2))


# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))







# ****************************************************** #
#### TRACK READS THROUGH THE PIPELINE #### 
# ****************************************************** #
# check number of reads that remained through each step in pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), 
                     rowSums(seqtab.nochim), rowSums(seqtab.nochim)/out[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "tabled", "nonchim", "%_of_start")
rownames(track) <- sample.names
head(track)

mean(track[,2] / track[,1], na.rm = T) # 0.642
mean(track[,3] / track[,2], na.rm = T) # 0.131
mean(track[,5] / track[,4], na.rm = T) # 0.979
mean(track[,6]) # 0.083



# *********************************** #


getN <- function(x) sum(getUniques(x))
track.2 <- cbind(out.2, sapply(mergers.2, getN), rowSums(seqtab.2), 
                 rowSums(seqtab.nochim.2), rowSums(seqtab.nochim.2)/out.2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track.2) <- c("input", "filtered", "merged", "tabled", "nonchim", "%_of_start")
rownames(track.2) <- sample.names
head(track.2)

mean(track.2[,2] / track.2[,1], na.rm = T) # 0.725
mean(track.2[,3] / track.2[,2], na.rm = T) # 0.915 ** The merging step was far more successful here
mean(track.2[,5] / track.2[,4], na.rm = T) # 0.849
mean(track.2[,6]) # 0.564 -- so this filtering retains way more than the first filtering... should be much more useful














# ****************************************************** #
#### GET OTU TABLE (genus level) #### 
# ****************************************************** #

# change NAs to unclassified
tax.fixNames <- tax
tax.fixNames[is.na(tax.fixNames)] <- "unclassified"
# combine/adjust genus names when appropriate
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Butyrivibrio_2"), "Genus"] <- "Butyrivibrio"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <- "Clostridium"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Corynebacterium_1"), "Genus"] <- "Corynebacterium"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <- "Escherichia"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Prevotella_2","Prevotella_6","Prevotella_7","Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Selenomonas_3","Selenomonas_4"), "Genus"] <- "Selenomonas"
tax.fixNames[ tax.fixNames[,"Genus"] %in% c("Spirochaeta_2"), "Genus"] <- "Spirochaeta"


# add other levels to unclassified genera to have unique values
tax.fixNames <- t(apply(tax.fixNames, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.fixNames), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                 paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                 all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.fixNames) <- colnames(tax)
# get otu table
otu.table <- seqtab.nochim
colnames(otu.table) <- unname(tax.fixNames[colnames(seqtab.nochim) , "Genus"])
non.bact <- tax.fixNames[ tax.fixNames[,"Kingdom"] %in% c(NA, "Eukaryota","Archaea"), "Genus"]
otu.table <- otu.table[ , ! colnames(otu.table) %in% non.bact ]
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table <- t(otu.table %*% sapply(unique(colnames(otu.table)),"==", colnames(otu.table)))
otu.table.rel <- apply(otu.table, 2, function(x) 100 * x/sum(x))
otu.table.rel[is.nan(otu.table.rel)] <- 0

write.csv(otu.table, "/users/tg/jwillis/SLL/Part_2/DADA2/otu_table.csv")
write.csv(otu.table.rel, "/users/tg/jwillis/SLL/Part_2/DADA2/otu_table_rel.csv")


# get tax table with rownames as genus values
tax.table <- unique(tax.fixNames)
rownames(tax.table) <- unname(tax.table[,"Genus"])
tax.table <- tax.table[ ! rownames(tax.table) %in% non.bact, ]
tax.table <- tax.table[ sort(rownames(tax.table)), ]
write.csv(tax.table, "/users/tg/jwillis/SLL/Part_2/DADA2/tax_table.csv")




# ****************************************************** #


# change NAs to unclassified
tax.fixNames.2 <- tax.2
tax.fixNames.2[is.na(tax.fixNames.2)] <- "unclassified"
# combine/adjust genus names when appropriate
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Butyrivibrio_2"), "Genus"] <- "Butyrivibrio"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <- "Clostridium"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Corynebacterium_1"), "Genus"] <- "Corynebacterium"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <- "Escherichia"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Prevotella_2","Prevotella_6","Prevotella_7","Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Selenomonas_3","Selenomonas_4"), "Genus"] <- "Selenomonas"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Spirochaeta_2"), "Genus"] <- "Spirochaeta"
tax.fixNames.2[ tax.fixNames.2[,"Genus"] %in% c("Treponema_2"), "Genus"] <- "Treponema"

# add other levels to unclassified genera to have unique values
tax.fixNames.2 <- t(apply(tax.fixNames.2, 1, function(x) {
  all.levels <- as.character(x)
  if ("unclassified" %in% all.levels) {
    new.row <- sapply(2:ncol(tax.fixNames.2), function(tl) ifelse(all.levels[tl]=="unclassified", 
                                                                paste(c("unclassified", all.levels[1:(tl-1)]), collapse = '.'),
                                                                all.levels[tl]))
    new.row <- c(all.levels[1], new.row)
  } else {
    all.levels
  }
  
}))
colnames(tax.fixNames.2) <- colnames(tax.2)
# get otu table
otu.table.2 <- seqtab.nochim.2
colnames(otu.table.2) <- unname(tax.fixNames.2[colnames(seqtab.nochim.2) , "Genus"])
non.bact.2 <- tax.fixNames.2[ tax.fixNames.2[,"Kingdom"] %in% c(NA, "Eukaryota","Archaea"), "Genus"]
otu.table.2 <- otu.table.2[ , ! colnames(otu.table.2) %in% non.bact.2 ]
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.table.2 <- t(otu.table.2 %*% sapply(unique(colnames(otu.table.2)),"==", colnames(otu.table.2)))
otu.table.rel.2 <- apply(otu.table.2, 2, function(x) 100 * x/sum(x))
otu.table.rel.2[is.nan(otu.table.rel.2)] <- 0

write.csv(otu.table.2, "/users/tg/jwillis/SLL/Part_2/DADA2/otu_table.2.csv")
write.csv(otu.table.rel.2, "/users/tg/jwillis/SLL/Part_2/DADA2/otu_table_rel.2.csv")


# get tax table with rownames as genus values
tax.table.2 <- unique(tax.fixNames.2)
rownames(tax.table.2) <- unname(tax.table.2[,"Genus"])
tax.table.2 <- tax.table.2[ ! rownames(tax.table.2) %in% non.bact.2, ]
tax.table.2 <- tax.table.2[ sort(rownames(tax.table.2)), ]
write.csv(tax.table.2, "/users/tg/jwillis/SLL/Part_2/DADA2/tax_table.2.csv")




















# ****************************************************** #
#### GET OTU TABLE (species level) #### 
# ****************************************************** #

# change NAs to unclassified
tax.fixNames.spec <- tax_spec
tax.fixNames.spec[is.na(tax.fixNames.spec)] <- "unclassified"
# combine/adjust genus names when appropriate
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Butyrivibrio_2"), "Genus"] <- "Butyrivibrio"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <- "Clostridium"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Corynebacterium_1"), "Genus"] <- "Corynebacterium"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <- "Escherichia"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Prevotella_2","Prevotella_6","Prevotella_7","Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames.spec[ tax.fixNames.spec[,"Genus"] %in% c("Selenomonas_3","Selenomonas_4"), "Genus"] <- "Selenomonas"

# make species value scientific name by combining with genus name if not unclassified
tax.fixNames.spec[,"Species"] <- sapply(rownames(tax.fixNames.spec), function(ro) {
  if (tax.fixNames.spec[ro,"Species"] == "unclassified") {
    tax.fixNames.spec[ro,"Species"]
  } else {
    paste(tax.fixNames.spec[ro,"Genus"], tax.fixNames.spec[ro,"Species"], sep = " ")
  }
})

# add other levels to unclassified genera to have unique values
tax.fixNames.spec <- t(apply(tax.fixNames.spec, 1, function(x) {
  all.levels.spec <- as.character(x)
  if ("unclassified" %in% all.levels.spec) {
    new.row.spec <- sapply(2:ncol(tax.fixNames.spec), function(tl) ifelse(all.levels.spec[tl]=="unclassified", 
                                                                          paste(c("unclassified", all.levels.spec[1:(tl-1)]), collapse = '.'),
                                                                          all.levels.spec[tl]))
    new.row.spec <- c(all.levels.spec[1], new.row.spec)
  } else {
    all.levels.spec
  }
  
}))
colnames(tax.fixNames.spec) <- colnames(tax_spec)
# get otu.spec table
otu.spec.table <- seqtab.nochim
colnames(otu.spec.table) <- unname(tax.fixNames.spec[colnames(seqtab.nochim) , "Species"])
non.bact.spec <- tax.fixNames.spec[ tax.fixNames.spec[,"Kingdom"] %in% c(NA, "Eukaryota","Archaea"), "Species"]
otu.spec.table <- otu.spec.table[ , ! colnames(otu.spec.table) %in% non.bact.spec]
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.spec.table <- t(otu.spec.table %*% sapply(unique(colnames(otu.spec.table)),"==", colnames(otu.spec.table)))
otu.spec.table.rel <- apply(otu.spec.table, 2, function(x) 100 * x/sum(x))
otu.spec.table.rel[is.nan(otu.spec.table.rel)] <- 0

write.csv(otu.spec.table, sprintf("%s/Part_2/DADA2/otu_spec_table.csv", home_dir))
write.csv(otu.spec.table.rel, sprintf("%s/Part_2/DADA2/otu_spec_table_rel.csv", home_dir))


# get tax table with rownames as species values
tax.spec.table <- unique(tax.fixNames.spec)
rownames(tax.spec.table) <- unname(tax.spec.table[,"Species"])
tax.spec.table <- tax.spec.table[ ! rownames(tax.spec.table) %in% non.bact.spec, ]
tax.spec.table <- tax.spec.table[ sort(rownames(tax.spec.table)), ]
write.csv(tax.spec.table, sprintf("%s/Part_2/DADA2/tax_spec_table.csv", home_dir))



# ****************************************************** #


# change NAs to unclassified
tax.fixNames.spec.2 <- tax_spec.2
tax.fixNames.spec.2[is.na(tax.fixNames.spec.2)] <- "unclassified"
# combine/adjust genus names when appropriate
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Butyrivibrio_2"), "Genus"] <- "Butyrivibrio"
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Clostridium_sensu_stricto_1"), "Genus"] <- "Clostridium"
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Corynebacterium_1"), "Genus"] <- "Corynebacterium"
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Escherichia/Shigella"), "Genus"] <- "Escherichia"
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Prevotella_2","Prevotella_6","Prevotella_7","Prevotella_9"), "Genus"] <- "Prevotella"
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Selenomonas_3","Selenomonas_4"), "Genus"] <- "Selenomonas"
tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Genus"] %in% c("Treponema_2"), "Genus"] <- "Treponema"

# make species value scientific name by combining with genus name if not unclassified
tax.fixNames.spec.2[,"Species"] <- sapply(rownames(tax.fixNames.spec.2), function(ro) {
  if (tax.fixNames.spec.2[ro,"Species"] == "unclassified") {
    tax.fixNames.spec.2[ro,"Species"]
  } else {
    paste(tax.fixNames.spec.2[ro,"Genus"], tax.fixNames.spec.2[ro,"Species"], sep = " ")
  }
})

# add other levels to unclassified genera to have unique values
tax.fixNames.spec.2 <- t(apply(tax.fixNames.spec.2, 1, function(x) {
  all.levels.spec <- as.character(x)
  if ("unclassified" %in% all.levels.spec) {
    new.row.spec <- sapply(2:ncol(tax.fixNames.spec.2), function(tl) ifelse(all.levels.spec[tl]=="unclassified", 
                                                                          paste(c("unclassified", all.levels.spec[1:(tl-1)]), collapse = '.'),
                                                                          all.levels.spec[tl]))
    new.row.spec <- c(all.levels.spec[1], new.row.spec)
  } else {
    all.levels.spec
  }
  
}))
colnames(tax.fixNames.spec.2) <- colnames(tax_spec.2)
# get otu.spec table
otu.spec.table.2 <- seqtab.nochim.2
colnames(otu.spec.table.2) <- unname(tax.fixNames.spec.2[colnames(seqtab.nochim.2) , "Species"])
non.bact.spec.2 <- tax.fixNames.spec.2[ tax.fixNames.spec.2[,"Kingdom"] %in% c(NA, "Eukaryota","Archaea"), "Species"]
otu.spec.table.2 <- otu.spec.table.2[ , ! colnames(otu.spec.table.2) %in% non.bact.spec.2 ]
# see this post for how to merge all the columns with the same column names:
#   https://stackoverflow.com/questions/11512441/combine-columns-in-matrix-having-same-column-name
otu.spec.table.2 <- t(otu.spec.table.2 %*% sapply(unique(colnames(otu.spec.table.2)),"==", colnames(otu.spec.table.2)))
otu.spec.table.rel.2 <- apply(otu.spec.table.2, 2, function(x) 100 * x/sum(x))
otu.spec.table.rel.2[is.nan(otu.spec.table.rel.2)] <- 0

write.csv(otu.spec.table.2, sprintf("%s/Part_2/DADA2/otu_spec_table.2.csv", home_dir))
write.csv(otu.spec.table.rel.2, sprintf("%s/Part_2/DADA2/otu_spec_table_rel.2.csv", home_dir))


# get tax table with rownames as species values
tax.spec.table.2 <- unique(tax.fixNames.spec.2)
rownames(tax.spec.table.2) <- unname(tax.spec.table.2[,"Species"])
tax.spec.table.2 <- tax.spec.table.2[ ! rownames(tax.spec.table.2) %in% non.bact.spec.2, ]
tax.spec.table.2 <- tax.spec.table.2[ sort(rownames(tax.spec.table.2)), ]
write.csv(tax.spec.table.2, sprintf("%s/Part_2/DADA2/tax_spec_table.2.csv", home_dir))

write.csv(tax.fixNames.spec.2, file = sprintf("%s/Part_2/DADA2/tax.fixNames.spec.2.csv", home_dir))







# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))

# *********** #








# ****************************************************** #
#### EVALUATE ACCURACY #### 
# ****************************************************** #
# In the mock community here, mixture of 20 known strains (ref seqs in HMP_MOCK.v35.fasta)
# Now compare seq variants produced here to expected composition of mock community

mocks <- c(sample.names[ grepl('HM', sample.names) ], 
           sample.names[ grepl('MZ', sample.names) ], 
           sample.names[ grepl('NTC', sample.names) ])
mock.ref <- getSequences(file.path(path.tut, "HMP_MOCK.v35.fasta"))
# mock.genera <- c("Acinetobacter","Actinomyces","Bacillus","Bacteroides","Clostridium_sensu_stricto_1","Deinococcus","Enterococcus",
#                  "Escherichia/Shigella","Helicobacter","Lactobacillus","Listeria","Neisseria","Propionibacterium","Pseudomonas",
#                  "Rhodobacter","Staphylococcus","Staphylococcus","Streptococcus","Streptococcus","Streptococcus")
# # changed Escherichia => Escherichia/Shigella and Clostridium => Clostridium_sensu_stricto_1
# # because that is how the names appear in the silva database and thus how they are assigned here
mock.genera <- c("Acinetobacter","Actinomyces","Bacillus","Bacteroides","Clostridium","Deinococcus","Enterococcus",
                 "Escherichia","Helicobacter","Lactobacillus","Listeria","Neisseria","Propionibacterium","Pseudomonas",
                 "Rhodobacter","Staphylococcus","Staphylococcus","Streptococcus","Streptococcus","Streptococcus")

# values found here: http://downloads.ihmpdcc.org/data/HMMC/HMPRP_sT1-Mock.pdf
mock.genera.copies <- c(10000,1000,100000,1000,100000,1000,1000,1000000,10000,10000,10000,10000,10000,100000,1000000,100000,1000000,
                        100000,1000000,1000)
mock.genera.copies <- sapply(mock.genera.copies, function(x) round(100 * x/sum(mock.genera.copies), 2))
mock.genera.even <- rep(10000, 20)
mock.genera.even <- sapply(mock.genera.even, function(x) round(100 * x/sum(mock.genera.even), 2))
mock.genera.copies <- t(cbind(mock.genera, mock.genera.even, mock.genera.copies))
colnames(mock.genera.copies) <- mock.genera.copies[1,]
mock.genera.copies <- mock.genera.copies[2:3,]
mock.genera.copies <- apply(mock.genera.copies, 2, function(x) as.numeric(x))
rownames(mock.genera.copies) = c("Even","Staggered")
mock.genera.copies <- t(mock.genera.copies %*% sapply(unique(colnames(mock.genera.copies)), "==", colnames(mock.genera.copies)))




unqs.mock.list <- list()
for (m in mocks) {
  unqs.mock <- otu.table[,m]
  m.unqs <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
  unqs.mock.list[[m]] <- m.unqs
  
  cat("DADA2 inferred", length(m.unqs), "sample sequences present in the Mock community:", m, ".\n")
  
  if (length(m.unqs) > 0) {
    # match.ref <- sum(sapply(names(m.unqs), function(x) any(grepl(x, mock.ref))))
    match.ref <- sum(names(m.unqs) %in% mock.genera)
    cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")
  } else{
    cat("\n")
  }
}

# Residual error rate here is 0%. Mothur found 34 OTUs...so far fewer FPs here



unqs.mock.list <- list()
for (m in mocks) {
  unqs.mock <- otu.table.2[,m]
  m.unqs <- sort(unqs.mock[unqs.mock > 0], decreasing=TRUE) # Drop ASVs absent in the Mock
  unqs.mock.list[[m]] <- m.unqs
  
  cat("DADA2 inferred", length(m.unqs), "sample sequences present in the Mock community:", m, ".\n")
  
  if (length(m.unqs) > 0) {
    # match.ref <- sum(sapply(names(m.unqs), function(x) any(grepl(x, mock.ref))))
    match.ref <- sum(names(m.unqs) %in% mock.genera)
    cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n\n")
  } else{
    cat("\n")
  }
}

# Residual error rate here is 0%. Mothur found 34 OTUs...so far fewer FPs here






# Assign taxonomy to mock samples
snm <- seqtab.nochim[mocks, ]
snm <- snm[ , colSums(snm) > 0]
tax.mocks <- assignTaxonomy(snm, sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)

# Add species
tax_spec.mocks <- addSpecies(tax.mocks, sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))




# Assign taxonomy to mock samples
snm.2 <- seqtab.nochim.2[mocks, ]
snm.2 <- snm.2[ , colSums(snm.2) > 0]
tax.mocks.2 <- assignTaxonomy(snm.2, sprintf("%s/silva_nr_v132_train_set.fa.gz",path.tut), multithread=TRUE)

# Add species
tax_spec.mocks.2 <- addSpecies(tax.mocks.2, sprintf("%s/silva_species_assignment_v132.fa.gz",path.tut))









# ********************************* #
# function to get scores for each type of mock community
get_mock_comp_score <- function(m, ot, otr) {
  
  # use sum of counts in given mock community as a weight for score
  mock.counts <- sum(ot[ , m ])
  # use relative abundance values to compare to known percentages in mock communities
  mock.rel <- otr[ , m ]
  # add known percentage for those genera that do not appear in a given sample
  not_present <- unique(mock.genera[ ! mock.genera %in% names(mock.rel) ])
  
  if (startsWith(m, "HM-782")) {
    # for the even mock community
    sco <- sum(sapply(names(mock.rel), function(x) {
      ifelse(x %in% mock.genera,
             # if one of the genera that should be present, get percent difference
             abs(mock.rel[x] - mock.genera.copies[x,"Even"]) / mock.genera.copies[x,"Even"],
             # else it is a false positive, just add its abundance
             mock.rel[x])
    }))
    sco <- sco + sum(mock.genera.copies[not_present,"Even"])
    
  } else if (startsWith(m, "HM-783")) {
    # for the staggered mock community
    sco <- sum(sapply(names(mock.rel), function(x) {
      ifelse(x %in% mock.genera,
             # if one of the genera that should be present, get percent difference
             abs(mock.rel[x] - mock.genera.copies[x,"Staggered"]) / mock.genera.copies[x,"Staggered"],
             # else it is a false positive, just add its abundance
             mock.rel[x])
    }))
    sco <- sco + sum(mock.genera.copies[not_present,"Staggered"])
    
  } else {
    # for the control communities
    sco <- sum(mock.rel)
  }
  
  return(sco * mock.counts)
}
# ********************************* #



total_counts <- sum( otu.table )
total_counts.2 <- sum( otu.table.2 )

# score is a comparison to the known values in the mock community samples
#   calculates percent difference for the 17 genera that should appear in these samples
#   adds full value of any mock genus that does not appear in the sample (as a penalty)
#   then also adds value for any genus that should not be in the sample (as a penalty)
# So the lower the score the better
score <- sum(sapply(mocks, get_mock_comp_score, otu.table, otu.table.rel )) / total_counts # 0.472112
score.2 <- sum(sapply(mocks, get_mock_comp_score, otu.table.2, otu.table.rel.2 )) / total_counts.2 # 0.07830165

# ****** so the second method is way better and these are the tables that I will use ******



# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))


















unclassified.db <- readRDS("/users/tg/jwillis/SLL/unclassified_identifier_db.rds")



# ****************************************************** #
#### UPDATE unclassified identifier database #### 
# ****************************************************** #
# the numbers in the identifier will be arbitrary, but will remain the same for each unique taxonomy to maintain
# consistency of naming, even between different datasets (OCD, SLL, etc.)
tls <- c("Phylum","Class","Order","Family","Genus","Species")
for (tl in tls) {
  # new table subset to be added to unclassified identifier database
  table.subset <- tax.spec.table.2[ startsWith(tax.spec.table.2[ , tl], "unclassified"), 1:(match(tl, tls)+1) ]
  
  # only update table if there are any unclassified taxa at given level
  if ( length(table.subset) > 0 ) {
    
    # combine the existing unclassified identifier database with new table subset
    # but if only 1 row is present in the subset, will appear as a character vector instead of matrix, must force matrix class
    if ( class(table.subset) == "character" ) {
      tmp_table <- rbind( unclassified.db[[ tl ]], t(as.matrix(table.subset)) )
    } else {
      tmp_table <- rbind( unclassified.db[[ tl ]], table.subset )
    }
    
    # reduce the expanded unclassified values at all levels
    tmp_table[ startsWith( tmp_table, "unclassified") ] <- "unclassified"
    
    # keep only those rows which are unique for the columns up to but excluding the current taxa level 
    # (which is how the unclassified values are being identified)
    if ( is.null( rownames(unique(tmp_table[ , 1:match(tl, tls)])) ) ) {
      # this case occurs for the Phylum level because it checks just unique values in Kingdom,
      # so instead of returning another table (as in the else below), 
      # it returns an unnamed character vector, so cannot identify rows
      tmp_table <- unique(tmp_table)
    } else {
      tmp_table <- tmp_table[rownames(unique(tmp_table[ , 1:match(tl, tls)])),]
    }
    
    # update rownames to numbers, 
    # will keep current order so that those that were already present will be labeled with the same number
    rownames(tmp_table) <- 1:nrow(tmp_table)
    # Use numbered rows plus letter of tax level to label the deepest level unclassified taxa
    tletter <- substr(tl, 1, 1)
    tmp_table[ , tl] <- sapply(rownames(tmp_table), function(x) sprintf("unclassified.%s%s",tletter,x))
    rownames(tmp_table) <- tmp_table[ , tl]
    
    # finally, update the table for the given taxonomic level in the unclassified identifier database
    unclassified.db[[ tl ]] <- tmp_table
    
  }
}
saveRDS(unclassified.db, "/users/tg/jwillis/SLL/unclassified_identifier_db.rds")











# ****************************************************** #
#### UPDATE OTU and TAX TABLES with new unclassified identifiers #### 
# ****************************************************** #

t0 <- Sys.time()

otus <- list()
otus_rel <- list()
taxTables <- list()


# ************************************ #
# function for converting taxa names to appropriate unclassified identifiers
get_unclass_IDs <- function(taxa, level) {
  # get taxa table up to the indicated level
  tmp_tax <- unique(tax.spec.table.2[ , 1:(match(level, tls)+1) ])
  
  uID <- sapply(taxa, 
                function(x) {
                  if (startsWith(tmp_tax[x, level], "unclassified")) {
                    # get vector of taxa up to one level above current based on the extended
                    taxonomy <- tmp_tax[ x, 1:match(level, tls) ]
                    # remove upper taxa levels from all unclassified names, 
                    # in order to match to the unclassified identifier database in next step
                    taxonomy <- sapply(taxonomy, function(y) strsplit(y, '\\.')[[1]][1])
                    # get correct unclassified identifier by checking where taxonomy matches a row in the table for the given level
                    new.unclass <- rownames( unclassified.db[[ level ]] )[ sapply(rownames(unclassified.db[[ level ]]), 
                                                                                  function(y) sum(unclassified.db[[ level ]][y, 1:match(level, tls)] == taxonomy) == length(taxonomy) ) ]
                    return( new.unclass )
                  } else {
                    return( tmp_tax[x, level] )
                  }
                })
  if (level != "Species" & length( uID[ uID %in% uID[duplicated(uID)] ] ) > 0) {
    # this convoluted bit must exist because in this and other data, there are 2 instances of 
    # families called "Family_XI", but from 2 different orders, causes problems down the line
    # so I must give further identifiers to any instance of taxa names that are duplicated at a given
    # level, but actually have different taxonomy at higher levels
    dups <- uID[uID %in% uID[duplicated(uID)]]
    lev.up.from.dups <- tmp_tax[ names(dups), ncol(tmp_tax)-1 ]
    uID[ uID %in% uID[duplicated(uID)] ] <- sprintf("%s.%s", dups, lev.up.from.dups)
  }
  return( unname(uID) )

}

# ************************************ #
# Start by including the Species tax and otu tables, update the rownames with unclassified identifiers
taxTables[["Species"]] <- tax.spec.table.2
rownames(taxTables[["Species"]]) <- get_unclass_IDs( rownames(taxTables[["Species"]]), "Species")
taxTables[["Species"]] <- unique(taxTables[["Species"]])
otus[["Species"]] <- otu.spec.table.2
rownames(otus[["Species"]]) <- get_unclass_IDs( rownames(otus[["Species"]]), "Species")
otus_rel[["Species"]] <- otu.spec.table.rel.2
rownames(otus_rel[["Species"]]) <- get_unclass_IDs( rownames(otus_rel[["Species"]]), "Species")








# ************************************ #
# then get otu_tables at each level by summing the values from same taxa at that level within the species otu table
for (tl in c("Genus","Family","Order","Class","Phylum")) {
  print(tl)
  
  ### update tax table at given level first
  taxTables[[ tl ]] <- unique(tax.spec.table.2[ , 1:(match(tl, tls)+1) ])
  rownames(taxTables[[ tl ]]) <- get_unclass_IDs( rownames(taxTables[[ tl ]]), tl)
  
  ### get counts of all taxa at given tl from Species counts table that are within each value of the given tl
  taxa <- taxTables[[ tl ]][ , tl]
  otus[[ tl ]] <- apply(otus[["Species"]], 2,
                        function(x) sapply(taxa,
                                           function(y) sum(x[ names( taxTables[["Species"]][,tl][ y == taxTables[["Species"]][,tl] ] )]
                                           )
                        )
  )
  # example just to verify by hand that values are correct (they are now:
  # otus[["Genus"]][ 101:105, 1:5 ]
  # sum(otus[["Species"]][,4][ names( taxTables[["Species"]][,"Genus"][ taxTables[["Species"]][,"Genus"]==taxTables[["Genus"]][,"Genus"]["F0058"] ] )])
  # otus[["Species"]][names( taxTables[["Species"]][,"Genus"][ taxTables[["Species"]][,"Genus"]==taxTables[["Genus"]][,"Genus"]["F0058"] ] ), 4]
  
  # get normalized values -- relative abundances
  otus_rel[[ tl ]] <- apply( otus[[ tl ]], 2, function(x) 100 * x/sum(x))

}






# ************************************ #
# Make species names proper scientific names (by making them <Genus species>)

rownames(otus[["Species"]]) <- unname(sapply(rownames(otus[["Species"]]), 
                                             function(x) {
                                               if (startsWith(x, "unclassified")) {
                                                 if (startsWith(taxTables[["Species"]][x,"Genus"], "unclassified")) {
                                                   paste("unclassified", x, sep = ' ')
                                                 } else {
                                                   paste(taxTables[["Species"]][x,"Genus"], x, sep = ' ')
                                                 }
                                               } else {
                                                 x
                                               }
                                               
                                             }))

rownames(otus_rel[["Species"]]) <- unname(sapply(rownames(otus_rel[["Species"]]), 
                                                 function(x) {
                                                   if (startsWith(x, "unclassified")) {
                                                     if (startsWith(taxTables[["Species"]][x,"Genus"], "unclassified")) {
                                                       paste("unclassified", x, sep = ' ')
                                                     } else {
                                                       paste(taxTables[["Species"]][x,"Genus"], x, sep = ' ')
                                                     }
                                                   } else {
                                                     x
                                                   }
                                                   
                                                 }))

rownames(taxTables[["Species"]]) <- unname(sapply(rownames(taxTables[["Species"]]), 
                                                  function(x) {
                                                    if (startsWith(x, "unclassified")) {
                                                      if (startsWith(taxTables[["Species"]][x,"Genus"], "unclassified")) {
                                                        paste("unclassified", x, sep = ' ')
                                                      } else {
                                                        paste(taxTables[["Species"]][x,"Genus"], x, sep = ' ')
                                                      }
                                                    } else {
                                                      x
                                                    }
                                                    
                                                  }))


# ************************************ #
# fix the problem with Family_XI, which is found in both the orders Bacillales and Clostridiales
for (tl in c("Species","Genus","Family")) {
  taxTables[[ tl ]][ , "Family" ] <- unname(sapply( 1:nrow(taxTables[[tl]]), function(x) {
    if (taxTables[[tl]][ x, "Family" ] == "Family_XI") {
      sprintf("%s.%s", taxTables[[tl]][ x, "Family" ], taxTables[[tl]][ x, "Order" ])
    } else {
      taxTables[[tl]][ x, "Family" ]
    }
  } ))
  
  if (tl == "Family" & "Family_XI" %in% rownames(taxTables[[ tl ]]))
    rownames(taxTables[[ tl ]]) <- unname(sapply( rownames(taxTables[[ tl ]]), function(x) {
      if (x == "Family_XI") {
        taxTables[[ tl ]][ x, "Family"]
      } else {
        x
      }
    } ))
}




# ************************************ #
# replace long unclassified names with correct unclassified IDs in each column at each level
for (tl in rev(tls)) {
  #  must go backwards starting with deepest levels because will rely on values within higher levels
  for (col_to_change in tls[ (match(tl, tls)):1 ]) {
    # for given column, check the rownames from the table of that tax level to see what the official label should be
    #   in the case that col_to_change is same as tl, will essentially just use that tables rownames
    taxTables[[ tl ]][ , col_to_change ] <- unname(sapply( unname(taxTables[[ tl ]][ , col_to_change ]), function(x)
      rownames(taxTables[[ col_to_change ]])[ taxTables[[ col_to_change ]][, col_to_change ] == x ] ))
  }
}

# ************************************ #




# save tax and otu table objects
saveRDS(taxTables, "/users/tg/jwillis/SLL/Part_2/DADA2/taxTables.rds")
saveRDS(otus, "/users/tg/jwillis/SLL/Part_2/DADA2/otus.rds")
saveRDS(otus_rel, "/users/tg/jwillis/SLL/Part_2/DADA2/otus_rel.rds")

print(Sys.time() - t0)
# ************************************ #





# save.image(file = sprintf("%s/DADA2_files/RData_files/dada2_SLL2.RData", home_dir))










# ****************************************************** #
#### CREATE PHYLOGENETIC TREE #### 
# ****************************************************** #

# From: https://f1000research.com/articles/5-1492/v2
# "Phylogenetic relatedness is commonly used to inform downstream analyses, especially the 
#  calculation of phylogeny-aware distances between microbial communities. The DADA2 sequence 
#  inference method is reference-free, so we must construct the phylogenetic tree relating the 
#  inferred sequence variants de novo. We begin by performing a multiple-alignment using the 
#  DECIPHER R package11."


# First randomly sampling 1 sequence from all sequences assigned to each species
#   so will end with vector the length of number of species
#   in theory, the sequences assigned to a given species should be similar enough
#     that they would align together in the total alignment (10526 sequences here)
#     but can now do this much more efficiently (only 489 now)
non.euk <- names(sapply(unique(tax.fixNames.spec.2[,"Species"]), 
                        function(x) x %in% rownames(otu.spec.table.2))[ sapply(unique(tax.fixNames.spec.2[,"Species"]), 
                                                                             function(x) x %in% rownames(otu.spec.table.2))])
species.seqs <- sapply(sort(non.euk), 
                       function(x) sample(rownames(as.data.frame(tax.fixNames.spec.2)[tax.fixNames.spec.2[,"Species"]==x,]),1))
# get proper ID for unclassified species
names(species.seqs) <- unname(unlist(get_unclass_IDs(names(species.seqs), "Species")))

# make proper scientific names by adding the Genus
# spec.only is a vector of species names, for which the name attributes are the full scientific names
spec.only <- sapply(rownames(taxTables[["Species"]]), function(x) strsplit(x,' ')[[1]][2])
names(species.seqs) <- unname(sapply(names(species.seqs), function(x) 
  if (startsWith(x, "unclassified")) {
    names(spec.only[spec.only==x])
  } else {
    x
  }))



library(DECIPHER)
alignment <- AlignSeqs(DNAStringSet(species.seqs), anchor=NA)



# "The phangorn R package is then used to construct a phylogenetic tree. 
#  Here we first construct a neighbor-joining tree, and then fit a GTR+G+I 
#  (Generalized time-reversible with Gamma rate variation) maximum likelihood 
#  tree using the neighbor-joining tree as a starting point."
library(phangorn)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

saveRDS(fitGTR, file = "/users/tg/jwillis/SLL/Part_2/R_objects/SLL2.fitGTR.rds")
















