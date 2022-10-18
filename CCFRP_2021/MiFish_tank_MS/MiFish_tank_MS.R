library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)

library(xlsx)
library(writexl)

path <- "~/Desktop/CCFRP_2021/MiFish_tank_MS/MiFish_tank_MS_fastq"
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# cutadapt
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
for(a in 1:9) {
b <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[a]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[a]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[a]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[a]]))
print(b)
}

cutadapt <- "/Users/tienly/.local/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, # we do not change the default allowed error rate of 0.1
                             "-m 1", # -m 1 discards reads having a length of zero bp after cutadapting to avoid problem
                             "-n 2", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}
# Count the presence of primers in the first cutdapt-ed sample to check if cutadapt worked as intended
for(a in 1:9) {
b <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[a]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[a]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[a]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[a]]))
print(b)
}

# dada2
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")
# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
# Set filtering parameter
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")
# Inspecting the returned dada-class object of the first sample:
dadaFs
dadaRs

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
# Inspect the merger data.frame of the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
# How many sequence variants were inferred?
dim(seqtab)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# create RDB
RDB <- read.xlsx("12S_marine_RDB_Gold.xlsx", 1)
write.table(RDB, file = "MiFish_RDB.txt", sep = ";",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# go to RDB_final-formatting.Rmd

taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/CCFRP_2021/MiFish_tank_MS/MiFish_RDB.txt", multithread=TRUE, tryRC=TRUE)
taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel
df <- as.data.frame(taxa)
df$ASV <- rownames(df)
write_xlsx(df, "~/Desktop/CCFRP_2021/MiFish_tank_MS/MiFish_tank_MS_taxonomy.xlsx")

# phyloseq
metadata <- read.csv("MiFish_MS_map.csv")
rownames(metadata) <- metadata$run_id
rownames(seqtab.nochim) <- metadata$run_id

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), tax_table(taxa))

plot_bar(ps, x="sample_ID", fill="Species") +
  # facet_grid(~protocol~site, scales="free") +
  facet_wrap(~protocol~site, scales="free") +
  labs(x = "Sample ID", y = "Abundance")  +
  labs(title = "MiFish Tank MS") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5))

# relative abundance
rownames(df) <- taxa_names(ps)
df$Count <- taxa_sums(ps)
df$Relative_abundance <- df$Count / sum(df$Count)
write_xlsx(df, "~/Desktop/CCFRP_2021/MiFish_tank_MS/MiFish_tank_MS_abun.xlsx")
