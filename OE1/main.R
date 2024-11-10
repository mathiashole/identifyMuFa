#!/usr/bin/env Rscript


# call all scripts
FILTER_SEQ <- "code/filter_keywords.sh"
GSCISSORS <- "code/Genomics_Scissors/gscissors.pl"
# SEQ_A <- "/home/usuario/Data_Rstudio/chop_genome/seq_attributes/seq_attributes.R"
# DISTRIBUTION <- "/home/usuario/Data_Rstudio/statistics_of_sequence/analyze_statistics.R"
#SEQ_A <- "/home/usuario/Data_Rstudio/statistics_of_sequence/render_quarto.R"

# Load necessary libraries
library(tidyverse)

# Create output directory if it does not exist
create_output_dir <- function(output_dir) {
  execution_path <- getwd()
  output_dir <- paste0(execution_path, "/", output_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  return(output_dir)
}

# Function to read input file and set column names
read_input <- function(input_file) {
  data <- read_tsv(input_file, col_names = FALSE)
  colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
  return(data)
}

transform_data <- function(data) {
  # Perform metadata transformations
  data <- mutate(data,
    gff_basename = basename(gff_file),
    no_gff_basename = str_remove(gff_basename, ".{4}$"),
    keyword_sum = paste(keyword1, keyword2, sep = "_"),
    filtred_name_gff = str_c("filtered_:", keyword_sum, ":_", gff_basename),
    non_filtred_name_gff = str_c("non-filtered_:", keyword2, ":_", gff_basename),
    out_gscissors = str_c("out_:", keyword_sum, ":_", no_gff_basename, ".fasta"),
    out_rest_gscissors = str_c("out_rest:", keyword2, ":_", no_gff_basename, ".fasta"),
    stat_fasta_feature = str_c("stat_", keyword_sum, "_", no_gff_basename, ".tsv"),
  )
  
  # Returns the transformed DataFrame
  return(data)
}

generate_commands <- function(data, output_dir) {
  # we use rowwise to apply the function to each row individually
  data <- rowwise(data) %>%
    mutate(
      filter_seq_command = str_c(FILTER_SEQ, " ", input_file, " ", output_dir),
      gscissors_command = str_c(GSCISSORS, " --fasta ", fasta_file, " --coordinates ", output_dir, "/", filtred_name_gff, " --format gff --output ", output_dir, "/", out_gscissors),
      gscissors_rest_command = str_c(GSCISSORS, " --fasta ", fasta_file, " --coordinates ", output_dir, "/", non_filtred_name_gff, " --format gff --output ", output_dir, "/", out_rest_gscissors)#,
      # fasta_feature_command = str_c(SEQ_A, " ", out_gscissors),
      # distribution_command = str_c(DISTRIBUTION, " ", stat_fasta_feature)
    )
  
  # Return DataFrame with generated commands
  return(data)
}

# Function to process each set of arguments
execution_module <- function(data, output_dir) {
      for (i in 1:nrow(data)) {
        cat("Processing FILTER_SEQ: ", data$filter_seq_command[i], "\n")
        system(data$filter_seq_command[i])
        # Check if FILTER_SEQ created the expected file
        path_file_fq <- paste0(output_dir,"/", data$filtred_name_gff[i])
        if (!file.exists(path_file_fq)) {
          cat("Error: FILTER_SEQ did not create the file", data$filtred_name_gff[i], "\n")
          next
        }

        # change directory from where you get the data!! DEBUGGING
        cat("Processing GSCISSORS: ", data$gscissors_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gs <- paste0(output_dir,"/", data$out_gscissors[i])
        system(data$gscissors_command[i])
        if (!file.exists(path_file_gs)) {
          cat("Error: GSCISSORS did not create the file", data$out_gscissors[i], "\n")
          next
        }

        cat("Processing GSCISSORS: ", data$gscissors_rest_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gs_rest <- paste0(output_dir,"/", data$out_rest_gscissors[i])
        system(data$gscissors_rest_command[i])
        if (!file.exists(path_file_gs_rest)) {
          cat("Error: GSCISSORS did not create the file", data$out_rest_gscissors[i], "\n")
          next
        }
        # # change directory from where you get the data!! DEBUGGING
        # cat("Processing SEQ_A: ", data$fasta_feature_command[i], "\n")
        # system(data$fasta_feature_command[i])
        # # change directory from where you get the data!! DEBUGGING
        # cat("Processing DISTRIBUTION: ", data$distribution_command[i], "\n")
        # system(data$distribution_command[i])

        # cat("Successfully processed", data$fasta_file[i], "and keywords:", data$keyword1[i], data$keyword2[i], "\n")
      }
}

# Check the number of arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: process_sequences.R <file_with_arguments.tsv>")
}

#####################
# Execution section #
#####################

input_file <- args[1]
output_dir <- "output_directory"  # Define your output directory
output_dir <- create_output_dir(output_dir)

data <- read_input(input_file)

# Apply the transformations
data_transformed <- transform_data(data)
path_data_transformed <- file.path(output_dir, paste0("data_transformed.tsv"))
write.table(data_transformed, file = path_genome_safe_tsv, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# Generate the commands
data_with_commands <- generate_commands(data_transformed, output_dir)
# Execution script
execution_module(data_with_commands, output_dir)