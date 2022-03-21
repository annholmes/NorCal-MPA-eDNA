library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)

# Getting ready #
path <- "~/Downloads/pilot"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME.1.fastq.gz and SAMPLENAME.2.fastq.gz
fnFs <- sort(list.files(path, pattern=".1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".2.fastq.gz", full.names = TRUE))

# Remove primers with cutadapt #
FWD <- "ACTGGGATTAGATACCCC" # forward primer sequence
REV <- "TAGAACAGGCTCCTCTAG" # reverse primer sequence

#Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input senquence
  require(Biostrings)
  dna <- DNAString(primer) # The Biostrings works with DNAString objects rather than character vector
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Calculate number of reads containing forward and reverse primer sequences (Only exact matches are found.)
# Only one set of paired end fastq.gz files will be checked (first sample in this case)
# This is sufficient, assuming all the files were created using the same library preparation.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
# Output interpretation
# FWD primer is found in the forward reads in its forward orientation
# FWD primer is found in the forward reads in its reverse-complement orientation (only 1)
# FWD primer is found in the reverse reads in its forward orientation (only 1)
# FWD primer is found in the reverse reads in its reverse-complement orientation


# Use cutadapt for primer removal:

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
# Trim FWD and its reverse-complement off of R1 (forward read)
R1.flags <- paste("-g", FWD, "-a", FWD.RC)
# Trim FWD and its reverse-complement off of R2 (reverse read)
R2.flags <- paste("-G", FWD.RC, "-A", FWD)
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, # we do not change the default allowed error rate of 0.1
                             "-n", 2, # -n 2 required to remove FWD and its reverse-complement from read
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
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
cutFs <- sort(list.files(path.cut, pattern = "1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "2.fastq.gz", full.names = TRUE))

# Check if forward and reverse files match:
if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")
# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

