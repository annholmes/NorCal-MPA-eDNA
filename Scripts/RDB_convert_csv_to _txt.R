#converts a csv with trimmed barcodes and taxonomy to text file format for use in dada2 pipelline
#taxonomy (K, P, C, O, F, G, S) must be correctly formatted in seven columns
#there is no column with both genus and species
#no metadata is present

#set working directory
setwd("/Users/annholmes/github/ann_holmes/CCFRP")

#import the csv file
RDB_seqs <-read.csv("12S_Marine_Seqs.csv", header = TRUE)

#remove any rows that do not have barcodes
RDB <- RDB_seqs[complete.cases(RDB_seqs), ]

#check the file 
> dim(RDB)

#export as text file
#file is named for the barcode and any specifics on the group of species (e.g. habitat), not the project
write.table(RDB, file = "12S_MiFish_Marine_RDB.txt", sep = ";",row.names = FALSE, col.names = FALSE, quote = FALSE)

