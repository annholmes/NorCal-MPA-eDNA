#NorCal MPA eDNA data analysis
#template for cutadapt + dada2 + cutadapt bioinformatics pipeline for multiple markers

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(xlsx)
library(writexl)
library(here)  #project-relative paths

#marker name and paths
marker_name <- "add_marker" #12S_MiFish, 12S_Riaz, 16S_Berry, CytB_MiSebastes
data_path <- here("data", marker_name)
results_path <- here("results", marker_name)
dir.create(results_path, recursive = TRUE, showWarnings = FALSE) #create folder if it doesn't exist already

#list fastq files for the marker of choice
fnFs <- sort(list.files(data_path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(data_path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
cat("Found", length(fnFs), "forward reads and", length(fnRs), "reverse reads for", marker_name, "\n")
if (length(fnFs) == 0 | length(fnRs) == 0) stop("No FASTQ files found in data_path.")

# cutadapt for <add marker name>
FWD <- "add_sequence" #forward primer sequence
REV <- "add sequence" #reverse primer sequence

#verify presence and orientation of these primers in the data
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

#calculate number of reads containing forward and reverse primer sequences (only exact matches are found)
#only one set of paired end fastq.gz files will be checked (first sample in this case)
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

for(a in 1:length(fnFs)) {
  b <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[a]]), 
             FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[a]]), 
             REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[a]]), 
             REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[a]]))
  print(b)
}

cutadapt <- "/opt/miniconda3/envs/cutadaptenv/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- here("data", marker_name, "cutadapt")  # cutadapt outputs inside the marker folder
if(!dir.exists(path.cut)) dir.create(path.cut, recursive = TRUE, showWarnings = FALSE)

filt_path <- file.path(path.cut, "filtered")
dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)

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

#count the presence of primers in the first cutadapt-ed sample to check if cutadapt worked as intended
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#dada2
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")
# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

#plotQualityProfile(cutFs[1:2])
#plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))  

#set marker-specific filtering parameters
#this marker is approx <XXX> bp
min_len <- <XXX> #add marker-specific min length (can also remove all together, if appropriate)
max_len <- <XXX> #add marker-specific max length (can also remove all together, if appropriate)
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = min_len, maxLen = max_len,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)

#check how many reads remain after filtering
out

set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")

#inspect the returned dada-class object of the first sample:
dadaFs
dadaRs

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
#inspect the merger data.frame of the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
#how many sequence variants were inferred?
dim(seqtab)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

rdb_file <- "rdb_name.txt"  #set manually for each marker

taxa <- assignTaxonomy(
  seqtab.nochim,
  here("resources", rdb_file),
  multithread = TRUE,
  tryRC = TRUE
)

taxa.print <- taxa  #remove sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel
df <- as.data.frame(taxa)
df$ASV <- rownames(df)
# relative abundance
df$Count <- colSums(seqtab.nochim)
df$Relative_abundance <- df$Count / sum(df$Count)
write_xlsx(df, file.path(results_path, paste0(marker_name, "_taxonomy.xlsx")))

saveRDS(seqtab, file = file.path(results_path, paste0(marker_name, "_ASV.rds")))  #ASV table

