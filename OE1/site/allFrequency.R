#!/usr/bin/env Rscript

# Call library
suppressPackageStartupMessages({
  library(seqinr)
  library(dplyr)
  library(parallel)

})

# ---- Advanced settings ----
generate_unique_filename <- function(type) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  paste0("COMBINED_", type, "_", timestamp, ".tsv")
}

# ---- Calculate function ----
# Function to count dinucleotides
count_dinucleotides <- function(sequence) {
  counts <- seqinr::count(sequence, 2)
  freqs <- counts / sum(counts)
  freqs <- round(freqs, 3)
  return(freqs)
}
# Function to count trinucleotides
count_trinucleotides <- function(sequence) {
  counts <- seqinr::count(sequence, 3)
  freqs <- counts / sum(counts)
  freqs <- round(freqs, 3)
  return(freqs)
}

nucleotide_frequency <- function(fasta_files, type_frequency = "dinucleotide") {
  results <- data.frame()  # Initialize the empty data frame for the results

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
  }) %>% 
    bind_rows() %>%
    replace(is.na(.), 0) #%>%
    # pivot_longer(cols = 3:22, names_to = "AminoAcid", values_to = "Freq") %>%
    # mutate(Freq = Freq * 100)
  
  return(results)
}

# ---- MAIN ----

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize variables
  aminoacid_files <- NULL
  dinucleotide_files <- NULL
  trinucleotide_files <- NULL
  output_base <- NULL
  generate_combined <- TRUE
  
  # Flag to know what we are reading
  current_mode <- NULL
  
  for (i in seq_along(args)) {
    if (args[i] %in% c("--aminoacid", "-aa")) {
      current_mode <- "aminoacid"
    } 
    else if (args[i] %in% c("--dinucleotide", "-di")) {
      current_mode <- "dinucleotide"
    }
    else if (args[i] %in% c("--trinucleotide", "-tri")) {
      current_mode <- "trinucleotide"
    }
    else if (args[i] %in% c("--output", "-o")) {
      output_base <- args[i + 1]
      i <- i + 1  # Skip value
    }
    else if (args[i] %in% c("--no-combined", "-nc")) {
      generate_combined <- FALSE
    }
    else {
      # Assign files according to the current mode
      if (!is.null(current_mode)) {
        if (current_mode == "aminoacid") {
          aminoacid_files <- c(aminoacid_files, args[i])
        } 
        else if (current_mode == "dinucleotide") {
          dinucleotide_files <- c(dinucleotide_files, args[i])
        }
        else if (current_mode == "trinucleotide") {
          trinucleotide_files <- c(trinucleotide_files, args[i])
        }
      }
    }
  }
  
    # Validate that only one type of analysis is specified
  analysis_types <- sum(!is.null(aminoacid_files), 
                     !is.null(dinucleotide_files), 
                     !is.null(trinucleotide_files))
  if (sum(analysis_types) != 1) {
    stop("Must specify exactly one analysis type: --aminoacid, --dinucleotide or --trinucleotide")
  }
  
  # Determine files and type of analysis
  if (!is.null(aminoacid_files)) {
    files <- aminoacid_files
    type <- "aminoacid"
  } else if (!is.null(dinucleotide_files)) {
    files <- dinucleotide_files
    type <- "dinucleotide"
  } else {
    files <- trinucleotide_files
    type <- "trinucleotide"
  }
  
  # Validate files
  existing_files <- files[file.exists(files)]
  if (length(existing_files) == 0) {
    stop("Error: None of the specified files exist")
  }
  
  return(list(
    type = type,
    files = existing_files,
    output_base = output_base,
    generate_combined = generate_combined
  ))
}


# args <- commandArgs(trailingOnly = TRUE)

# if (length(args) < 2) {
#   stop("Usage: ./calculate_nucleotide_frequencies.R <dinucleotide|trinucleotide> <fasta_file1> [fasta_file2 ...]")
# }

# type <- args[1]
# files <- args[-1]

# result <- nucleotide_frequncy(files, type_frequency = type)

# outfile <- paste0(type, "_frequencies.tsv")
# write.table(result, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
# cat("Output written to", outfile, "\n")
