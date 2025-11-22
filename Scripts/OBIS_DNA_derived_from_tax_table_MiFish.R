#eDNA to OBIS pipeline in R 
#converts taxa to WoRMS accepted names and uses robis for validation

library(tidyverse)
library(here)
library(digest) #for MD5 hashes
library(worrms) #WoRMS taxonomic alignment
library(robis) #validated tables for OBIS
library(dplyr)
library(purrr)
library(stringr)

#marker and paths
dataset_name <- "NorCal_MPA_eDNA"
marker_name <- "12S_MiFish"
results_path <- here("results", marker_name)
OBIS_path <- here("OBIS", marker_name)
dir.create(OBIS_path, recursive = TRUE, showWarnings = FALSE)

#marker and paths: automatically find files
tax_file  <- list.files(results_path, pattern = paste0(marker_name, "_taxonomy.*\\.csv$"), full.names = TRUE)
asv_file  <- list.files(results_path, pattern = paste0(marker_name, "_ASV_table.*\\.csv$"), full.names = TRUE)
meta_file <- list.files(results_path, pattern = paste0("NorCal_MPA_eDNA_metadata_cores_", marker_name, ".*\\.csv$"), full.names = TRUE)

# Read tables
tax_table  <- read_csv(tax_file)
asv_table  <- read_csv(asv_file)
meta_table <- read_csv(meta_file)

#add lowest rank column
# Define taxonomic hierarchy (lowest → highest)
rank_cols <- c("scientificName", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")

rank_names <- c(scientificName = "species",
                Genus          = "genus",
                Family         = "family",
                Order          = "order",
                Class          = "class",
                Phylum         = "phylum",
                Kingdom        = "kingdom")

get_taxon_rank <- function(row) {
  sn <- row[["scientificName"]] #check scientificName first
  if (!is.na(sn) && sn != "") {
    #two words → species
    if (grepl("^\\S+\\s+\\S+$", sn)) {
      return("species")
    }
    # one word → genus
    if (grepl("^\\S+$", sn)) {
      return("genus")
    }
  }
  #fallback to next available taxonomic column
  for (col in rank_cols[-1]) {  # skip scientificName
    if (!is.na(row[[col]]) && row[[col]] != "") {
      return(rank_names[[col]])
    }
  }
}

#apply to tax table
tax_table$taxonRank <- apply(tax_table, 1, get_taxon_rank)

#fill in missing taxonomy in tax table (i.e., replace NAs with lowest rank availale)
#columns to fill (adjust to your actual taxonomy columns)
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","scientificName")

#convert empty strings to NA if needed
tax_table <- tax_table %>%
  mutate(across(all_of(tax_cols), ~na_if(.x, "")))

#function to fill NA to the right within a row
fill_row_right <- function(x) {
  x <- as.character(x)
  last <- NA_character_
  for (i in seq_along(x)) {
    if (!is.na(x[i])) {
      last <- x[i]
    } else if (is.na(x[i]) && !is.na(last)) {
      x[i] <- last
    }
    # if both are NA, leave as NA
  }
  x
}

#apply row-wise to the selected columns
tax_table[tax_cols] <- t(apply(tax_table[tax_cols], 1, fill_row_right))

#fill verbatimScientificName
tax_table <- tax_table %>%
  mutate(
    verbatimScientificName = coalesce(verbatimScientificName, scientificName)
  )

#rename existing scientificName
tax_table <- tax_table %>%
  rename(scientificName_orig = scientificName)

#only look for each name once
unique_names <- unique(tax_table$scientificName_orig)

#function for WoRMS lookup
get_worms_exact <- function(name) {
  rec <- try(wm_records_name(name), silent = TRUE)
  
  if (inherits(rec, "try-error") || length(rec) == 0) {
    return(tibble(
      queried_name = name,
      AphiaID = NA_integer_,
      valid_name = NA_character_,
      scientificNameAuthorship = NA_character_   # <- add this
    ))
  }
  
  rec_df <- try(as.data.frame(rec), silent = TRUE)
  if (inherits(rec_df, "try-error")) {
    return(tibble(
      queried_name = name,
      AphiaID = NA_integer_,
      valid_name = NA_character_,
      scientificNameAuthorship = NA_character_   # <- add this
    ))
  }
  
  exact <- rec_df %>% filter(scientificname == name)
  
  if (nrow(exact) == 0) {
    return(tibble(
      queried_name = name,
      AphiaID = NA_integer_,
      valid_name = NA_character_,
      scientificNameAuthorship = NA_character_   # <- add this
    ))
  }
  
  row <- exact[1, ]
  
  tibble(
    queried_name = name,
    AphiaID = row$AphiaID,
    valid_name = row$valid_name,
    scientificNameAuthorship = row$authority      # <- rename authority
  )
}

#run lookup
unique_names <- unique(tax_table$scientificName_orig)

worms_lookup <- map_dfr(unique_names, get_worms_exact)

#add MD5 hash as sequenceID to tax table
tax_table <- tax_table %>%
  mutate(
    DNA_sequence = DNA_sequence,  # keep original column
    sequenceID = sapply(DNA_sequence, digest, algo = "md5")
  ) %>%
  relocate(sequenceID, .before = 1)  #puts sequenceID into column 1

#join WoRMS results back to tax_table
tax_table <- tax_table %>%
  left_join(worms_lookup, by = c("scientificName_orig" = "queried_name")) %>%
  mutate(
    scientificName = valid_name
  ) %>%
  select(
    sequenceID,
    scientificName,
    AphiaID,
    scientificNameAuthorship,
    Kingdom, Phylum, Class, Order, Family, Genus,
    verbatimScientificName,
    taxonRank,
    DNA_sequence
  )

#define non-marine taxa to flag
#likely detected from trace contamination
non_marine <- c("Cyprinus carpio", #only non-marine taxon with AphiaID in WoRMS
                "Homo sapiens", "Sus scrofa",
                "Capra hircus", "Meleagris gallopavo", "Felis catus",
                "Gallus gallus", "Canis familiaris")

#fill scientificName from verbatimScientificName where it's NA and add non-marine flag
tax_table <- tax_table %>%
  mutate(
    scientificName = ifelse(is.na(scientificName), verbatimScientificName, scientificName),
    nonMarineFlag = scientificName %in% non_marine
  )%>%
  relocate(nonMarineFlag, .after = scientificName)

#exclude from OBIS
tax_table <- tax_table %>%
  mutate(
    excludeFromOBIS = nonMarineFlag |
      (Class == "Mammalia" & taxonRank == "class") |
      (Class == "Actinopteri" & taxonRank == "class") |
      (Class == "Polychaeta" & taxonRank == "class") |
      (Class == "Mollusca" & taxonRank == "class")
  )

tax_table <- tax_table %>%
  mutate(
    sequence_remarks = ifelse(
      excludeFromOBIS,
      "Excluded from OBIS due to suspected trace contamination or insufficient taxonomic resolution",
      NA_character_
    )
  )

#create final DNA-derived OBIS table
DNA_derived_MiFish <- tax_table %>%
  mutate(
    env_broad_scale        = "marine biome [ENVO:00000447]",
    env_local_scale        = "coastal water [ENVO:00001250]",
    env_medium             = "waterborne particulate matter [ENVO:01000436]",
    nucl_acid_ext          = "Dneasy PowerWater Kit",
    nucl_acid_amp          = "2-step PCR",
    target_gene            = "12S rRNA",
    pcr_primer_forward     = "GTCGGTAAAACTCGTGCCAGC",
    pcr_primer_reverse     = "CATAGTGGGGTATCTAATCCCAGTTTG",
    pcr_primer_name_forward = "MiFish-U-F",
    pcr_primer_name_reverse = "MiFish-U-R",
    pcr_primer_reference   = "Miya et al. 2015 https://doi.org/10.1098/rsos.150088",
    seq_meth               = "Illumina MiSeq 2x150",
    lib_layout             = "paired",
    otu_class_appr         = "dada2; 1.36.0; ASV"
  ) %>%
  select(
    sequenceID, scientificName, AphiaID, scientificNameAuthorship,
    Kingdom, Phylum, Class, Order, Family, Genus,
    verbatimScientificName, taxonRank, DNA_sequence,
    env_broad_scale, env_local_scale, env_medium,
    nucl_acid_ext, nucl_acid_amp, target_gene,
    pcr_primer_forward, pcr_primer_reverse,
    pcr_primer_name_forward, pcr_primer_name_reverse,
    pcr_primer_reference, seq_meth, lib_layout, otu_class_appr,
    everything()
  )

#save the final file to csv
write_csv(DNA_derived_MiFish, here("OBIS", "12S_MiFish", "OBIS_DNA_derived_NorCal_MPA_eDNA_MiFish.csv"))

