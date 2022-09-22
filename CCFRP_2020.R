library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)

# Getting ready #
path <- "~/Downloads/pilot"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME.1.fastq.gz and SAMPLENAME.2.fastq.gz
fnFs <- sort(list.files(path, pattern = ".1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".2.fastq.gz", full.names = TRUE))

# Remove primers with cutadapt #
FWD <- "GTCGGTAAAACTCGTGCCAGC" # forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # reverse primer sequence
#Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Calculate number of reads containing forward and reverse primer sequences (Only exact matches are found.)
# Only one set of paired end fastq.gz files will be checked (first sample in this case)
# This is sufficient, assuming all the files were created using the same library preparation.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
# Output interpretation
# FWD primer is found in the forward reads in its forward orientation
# FWD primer is found in the reverse reads in its reverse-complement orientation 
# REV primer is found in the reverse reads in its forward orientation 
# REV primer is found in the forward reads in its reverse-complement orientation

# Use cutadapt for primer removal: #

# This requires prior installation of cutadapt on your machine via python, anaconda, etc.
# You can do this via python's pip function:
# Install python on your machine
# Run "pip install cutadapt" in command line

# Define the path to cutadapt.exe file:
cutadapt <- "/Library/Frameworks/Python.framework/Versions/3.9/bin/cutadapt"
# To see if this worked and cutadapt has indeed been found, check installed version of cutadapt:
system2(cutadapt, args = "--version") # Run shell commands from R. See if R recognizes cutadapt and shows its version
 
# Create output filenames for the cupadapt-ed files
# Define the parameters for the cutadapt command
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
 
# Run cutadapt
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, # we do not change the default allowed error rate of 0.1
                             "-m 1", # -m 1 discards reads having a length of zero bp after cutadapting to avoid problem
                             "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
# Count the presence of primers in the first cutdapt-ed sample to check if cutadapt worked as intended
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

## The primer-devoid sequence files are now ready to be analyzed. ##

# Similar to the earlier steps of reading in FASTQ files, read in the names of the cutadapt-ed FASTQ files.
# Get the matched lists of forward and reverse fastq files.

# Forward and reverse fastq filenames have the format (make sure file name patterns are correct):
cutFs <- sort(list.files(path.cut, pattern = ".1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = ".2.fastq.gz", full.names = TRUE))

# Check if forward and reverse files match:
if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")
# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

# Generate quality profile plots for our reads
# In case we have more than 20 fastq files, the following command will randomly choose 20 files to be plotted
if(length(cutFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(cutFs) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) #grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
}
fwd_qual_plots
rev_qual_plots

# Filter and trim #
# Assigning the directory for the filtered reads to be stored in
filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))
# Set filtering parameter
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

# Error model generation #
# dada2 learns the specific error-signature of our dataset
# This is why files from separate sequencing runs have to be processed separately till later on
set.seed(100) # set seed to ensure that randomized steps are replicable
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates:
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Apply the dada2's core sequence-variant inference algorithm: #
# Set pool = TRUE to allow information to be shared across samples
# This makes it easier to resolve rare variants which occur just once or twice in one sample but a few more times across samples
# This will increase computational time (most likely problematic once large data sets are analyzed)
# An alternative is pseudo-pooling, an intermediate solution with medium computational time
dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")

# Inspecting the returned dada-class object of the first sample:
dadaFs
dadaRs

# Merge the forward and reverse reads together to obtain the full denoised sequence #
# Adjust the minimum overlap (default = 12) and maximum mismatch allowed (e.g = 1) if necessary
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 10, maxMismatch = 1, verbose = TRUE)

# Inspect the merger data.frame of the first sample
head(mergers[[1]])

# Construct an amplicon sequence variant table (ASV) table #
seqtab <- makeSequenceTable(mergers)
# How many sequence variants were inferred?
dim(seqtab)

# Remove chimeras #
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
# Frequency of chimeric sequences
sum(seqtab.nochim)/sum(seqtab)

# Create a table to track read numbers throughout the pipeline: #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# Assign taxonomy #
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/12S_assignTaxonomy_formatted.fa", multithread=TRUE)

# Inspect the assignment result:
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
taxa.print
