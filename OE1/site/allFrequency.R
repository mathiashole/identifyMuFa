#!/usr/bin/env Rscript

# Call library
suppressPackageStartupMessages({
  library(seqinr)
  library(dplyr)
  library(parallel)

})

# Funciones

# ---- Advanced settings ----
generate_unique_filename <- function(type) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  paste0("COMBINED_", type, "_", timestamp, ".tsv")
}

# Function to count dinucleotides
count_dinucleotides <- function(sequence) {
  counts <- seqinr::count(sequence, 2)
  freqs <- counts / sum(counts)
  freqs <- round(freqs, 3)
  return(freqs)
}

# ---- Calculate function ----

# Function to count trinucleotides
count_trinucleotides <- function(sequence) {
  counts <- seqinr::count(sequence, 3)
  freqs <- counts / sum(counts)
  freqs <- round(freqs, 3)
  return(freqs)
}

nucleotide_frequncy <- function(fasta_files, type_frequency = "dinucleotide") {
  sequences <- read.fasta(fasta_file)
  base_name <- tools::file_path_sans_ext(basename(fasta_file))
  
  results <- data.frame()

  # Iterate over each FASTA file in the folder
  for (fasta_file in fasta_files) {
    sequences <- read.fasta(fasta_file)
    ids <- names(sequences)
    
    # Base name of file
    base_name <- basename(fasta_file)
    # Regular expresion for extract diferents parts
    extracted_names <- sub(".*_([^_]+)\\.fasta$", "\\1", base_name)
    
    # Get dinucleotide or trinucleotide frequencies for each sequence within a file
    for (i in seq_along(sequences)) {
      seq <- getSequence(sequences[[i]])

      if (type_frequency == "trinucleotide") {
        freqs <- count_trinucleotides(seq)
      } else if (type_frequency == "dinucleotide") {
        freqs <- count_dinucleotides(seq)
      }
      
      # Create a row with file name, ID and frequencies
      row <- c(File = extracted_names, ID = ids[i], freqs)
      # row <- data.frame(File = extracted_names, ID = ids[i], t(freqs), stringsAsFactors = FALSE)

      results <- results %>% bind_rows(row)
    }
  }

  return(results)
}

# Your exact amino acid frequency function (preserved)
calculate_amino_acid_frequencies <- function(file) {
  fasta <- read.fasta(file, seqtype = "AA")
  
  amino_acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  results <- lapply(names(fasta), function(seq_id) {
    seq <- toupper(getSequence(fasta[[seq_id]], as.string = TRUE)[[1]])
    aa_counts <- table(factor(strsplit(seq, "")[[1]], levels = amino_acids))
    total_aa <- sum(aa_counts)
    
    freq_df <- data.frame(AminoAcid = names(aa_counts), 
                         Frequency = as.vector(aa_counts) / total_aa)
    
    freq_wide <- freq_df %>%
      pivot_wider(names_from = AminoAcid, values_from = Frequency,
                 values_fill = list(Frequency = 0)) %>%
      mutate(File = basename(file), Sequence_ID = seq_id) %>%
      relocate(File, Sequence_ID)
    
    return(freq_wide)
  })
  
  combined <- bind_rows(results) %>%
    replace(is.na(.), 0)
  
  # Apply your filename transformation
  combined$File <- gsub(".*_([A-Za-z0-9]+)_translated\\.fasta", "\\1", combined$File)

  combined_long <- combined %>%
  pivot_longer(cols = c(3:22),  # Especifica las columnas de aminoácidos
               names_to = "AminoAcid",  # Nombre de la nueva columna para los nombres de los aminoácidos
               values_to = "Freq") %>%  # Nombre de la nueva columna para las frecuencias
  mutate(Freq = Freq * 100)  # Convertir frecuencias a porcentaje
  
  return(combined_long)
}

# ---- MAIN ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: ./calculate_nucleotide_frequencies.R <dinucleotide|trinucleotide> <fasta_file1> [fasta_file2 ...]")
}

type <- args[1]
files <- args[-1]

result <- nucleotide_frequncy(files, type_frequency = type)

outfile <- paste0(type, "_frequencies.tsv")
write.table(result, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Output written to", outfile, "\n")
