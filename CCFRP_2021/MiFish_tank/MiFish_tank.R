library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(psadd)

library(xlsx)
library(writexl)


# Getting ready #
path <- "~/Desktop/CCFRP_2021/MiFish_tank/MiFish_tank_fastq"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Remove primers with cutadapt #
FWD <- "GTCGGTAAAACTCGTGCCAGC" # forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # reverse primer sequence

#Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Calculate number of reads containing forward and reverse primer sequences (Only exact matches are found.)
# Only one set of paired end fastq.gz files will be checked (first sample in this case)
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# Use cutadapt for primer removal
# Define the path to cutadapt.exe file
cutadapt <- "/Users/tienly/.local/bin/cutadapt"
system2(cutadapt, args = "--version")

# Create output filenames for the cupadapt-ed files
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Define the parameters for the cutadapt command
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# R1.flags <- paste("-g", FWD, "-a", REV.RC)
# R2.flags <- paste("-G", REV, "-A", FWD.RC)
R1.flags <- paste("-g", FWD, "-g", REV, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-G", FWD, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, # we do not change the default allowed error rate of 0.1
                             "-m 1", # -m 1 discards reads having a length of zero bp after cutadapting to avoid problem
                             "-n 2", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}
# Count the presence of primers in the first cutdapt-ed sample to check if cutadapt worked as intended
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
# The primer-devoid sequence files are now ready to be analyzed

# dada2 #
# Similar to the earlier steps of reading in FASTQ files, read in the names of the cutadapt-ed FASTQ files.
# Get the matched lists of forward and reverse fastq files.
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")

# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

# Generate quality profile plots for our reads
# plotQualityProfile(cutFs[1:2])
# plotQualityProfile(cutRs[1:2])

# Filter and trim 
# Assigning the directory for the filtered reads to be stored in
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
# Set filtering parameter
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

# Error model generation
set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates:
# plotErrors(errF, nominalQ = TRUE)
# plotErrors(errR, nominalQ = TRUE)

# Apply the dada2's core sequence-variant inference algorithm
dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")
# Inspecting the returned dada-class object of the first sample:
dadaFs
dadaRs

# Merge the forward and reverse reads together to obtain the full denoised sequence
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
# Inspect the merger data.frame of the first sample
head(mergers[[1]])

# Construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

# Create a table to track read numbers throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# create RDB if needed
# RDB <- read.xlsx("12S_marine_RDB_Gold.xlsx", 1)
# write.table(RDB, file = "MiFish_RDB.txt", sep = ";",
#            row.names = FALSE, col.names = FALSE, quote = FALSE)
# go to RDB_final-formatting.Rmd

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/CCFRP_2021/MiFish_tank/MiFish_RDB.txt", multithread=TRUE, tryRC=TRUE)
# Inspect the assignment result
taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel #
df <- as.data.frame(taxa)
df$ASV <- rownames(df)
# relative abundance
df$Count <- colSums(seqtab.nochim)
df$Relative_abundance <- df$Count / sum(df$Count)
write_xlsx(df, "~/Desktop/CCFRP_2021/MiFish_tank/MiFish_tank_tax.xlsx")

# phylogenetic tree #
library(DECIPHER)
library(phangorn)
theme_set(theme_bw())

# Extract sequences from DADA2 output
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences

# Build a neighbor-joining tree
# then fit a maximum likelihood tree
# using the neighbor-joining tree as a starting point

# Run sequence alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
#Change sequence alignment output into a phyDat structure
phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
# Create distance matrix
dm <- dist.ml(phang_align)
# Perform neighbor joining
treeNJ <- NJ(dm)  # note, tip order != sequence order
# Internal maximum likelihood
fit = pml(treeNJ, data=phang_align)
# negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)

# phyloseq #
metadata <- read.csv("MiFish_tank_map.csv")
rownames(metadata) <- metadata$run_id
rownames(seqtab.nochim) <- metadata$run_id

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), tax_table(taxa), phy_tree(fitGTR$tree))
# Filter by length
ps.filt <- prune_taxa(nchar(taxa_names(ps)) >= 100 & nchar(taxa_names(ps)) <= 200, ps)
# Filter by number of reads
ps.filt <- prune_taxa(taxa_sums(ps.filt) >= 100, ps.filt) 

# Plot
plot_tree(ps.filt, color = "site", shape = "site", label.tips = "Species")

plot_bar(ps.filt, x="sample_ID", fill="Species") +
  # facet_grid(~protocol~site, scales="free") +
  facet_wrap(~protocol~site, scales="free") +
  labs(x = "Sample ID", y = "Abundance")  +
  labs(title = "MiFish Tank (filtered)") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5))

# Sample type
ps.nc <- subset_samples_no_zero(ps.filt, site == "pcr")
plot_tree(ps.nc, color="site", shape = "site", label.tips = "Species",
          ladderize = TRUE, plot.margin=0.1)

ps.filt1 <- subset_samples_no_zero(ps.filt, site != "pcr")
ps.MC <- subset_samples_no_zero(ps.filt1, protocol == "Cal eDNA")
ps.MS <- subset_samples_no_zero(ps.filt1, protocol == "SuperFi")

df <- as.data.frame(tax_table(ps.nc))
df$ASV <- rownames(df)
write_xlsx(df, "~/Desktop/CCFRP_2021/MiFish_tank/nc.xlsx")

df <- as.data.frame(tax_table(ps.MC))
df$ASV <- rownames(df)
write_xlsx(df, "~/Desktop/CCFRP_2021/MiFish_tank/mc.xlsx")

df <- as.data.frame(tax_table(ps.MS))
df$ASV <- rownames(df)
write_xlsx(df, "~/Desktop/CCFRP_2021/MiFish_tank/ms.xlsx")

a <- as.data.frame(taxa_sums(ps.filt))
a$ASV <- rownames(a)
