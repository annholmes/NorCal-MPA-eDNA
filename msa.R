#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msa")
library("msa")
library(microseq)
library(writexl)

# read in sequence file 
# all sequences must be in a single fasta file
seqs <- Biostrings::readDNAStringSet("sebastes_cytb/sebastes_cytb.fasta")
seqs

# create alignment using muscle algorithm
alignment <- msa(seqs, method="Muscle")
#print alignment results
print(alignment, show="complete")

# create fasta file of alignment
writeXStringSet(unmasked(alignment), file="sebastes_msa.fasta")

# trim multiple sequence alignment
# gap.end specifies the fraction of gaps tolerated at the end of the alignment
alignment <- readFasta("sebastes_msa.fasta")
print(str_length(alignment$Sequence))
msa_trimmed <- msaTrim(alignment, gap.end = 0, gap.mid = 0)
print(str_length(msa_trimmed$Sequence))

# export trimmed sequences as an excel file
write_xlsx(msa_trimmed, "~/Desktop/msa_trimmed.xlsx")
