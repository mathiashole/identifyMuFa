#!/usr/bin/env Rscript

# -------------------------------------------
# EXTERNAL MODULES
# -------------------------------------------

# External scripts called during the flow
FILTER_SEQ        <- "code/filter_keywords.sh"
GSCISSORS         <- "code/Genomics_Scissors/gscissors.pl"
OVERLAPPINGSHAIVE <- "code/overlappingShaive.R"
BOTHBLAST         <- "code/bothBlast.sh"
SPDIFFSIZE        <- "code/spDiffSize.R"
GORF              <- "code/gORF.sh"
ALLBLAST          <- "code/allBlast.sh"
BREFINER          <- "code/bRefiner.sh"
CLASSIFIER        <- "code/classifier.R"
# MEANSEQ        <- "code/meanSeq.sh"

# -------------------------------------------
# AUXILIARY FUNCTIONS
# -------------------------------------------

# Create output directory if it does not exist
create_output_dir <- function(output_dir) {
  execution_path <- getwd()
  output_dir <- paste0(execution_path, "/", output_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  return(output_dir)
}

read_input <- function(input_file) {
  read.delim(input_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
}

# Function to read input file and set column names
read_input <- function(input_file) {
  # data <- read_tsv(input_file, col_names = FALSE)
  # data <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE)
  data <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

  # Determine the number of columns
  num_cols <- ncol(data)
  # num_cols <- max(colSums(data != "" & !is.na(data))) 

  # Assign names based on the number of columns
#   if (num_cols == 6){
#     colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length", "sequence_file")
#   } else if (num_cols == 5){
#     colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length")
#   } else if (num_cols == 4) {
#     colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
#   } else if (num_cols == 3) {
#     colnames(data) <- c("fasta_file", "sequence_file", "length")
#   } else if (num_cols ==2) {
#     colnames(data) <- c("fasta_file", "sequence_file")
#   } else {
#     stop("Unexpected number of columns in input file.")
#   }

#   # return(list(data = data, num_cols = num_cols))
#   return(data)
# }

assign_column_names <- function(data, mode) {
  n <- ncol(data)

  if (mode == "hybrid") {
    if (n == 6) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length", "sequence_file")
    } else if (n == 5) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "sequence_file")
    } else {
      stop("Unexpected format in HYBRID mode")
    }
  } else if (mode == "gff") {
    if (n == 5) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length")
    } else if (n == 4) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
    } else {
      stop("Unexpected format in GFF mode")
    }
  } else if (mode == "no_gff") {
    if (n == 3) {
      colnames(data) <- c("fasta_file", "sequence_file", "length")
    } else if (n == 2) {
      colnames(data) <- c("fasta_file", "sequence_file")
    } else {
      stop("Unexpected format in NO_GFF mode")
    }
  }

  return(data)
}

# -------------------------------------------
# PARSEO DE ARGUMENTOS MANUAL
# -------------------------------------------

# Check the number of arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: process_sequences.R <file_with_arguments.tsv>")
}

# -------------------------------------------
# MAIN EXECUTION
# -------------------------------------------

input_file <- args[1]
# output_dir <- "output_directory"  # Define your output directory
# output_dir <- create_output_dir(output_dir)
# dir.create(file.path(output_dir, "blast_result"), showWarnings = FALSE, recursive = TRUE)

data <- read_input(input_file)
num_cols <- ncol(data)

if (num_cols == 3 || num_cols == 2) {
  
  source("code/moduleNoGff.R")

  print("without gff")
  output_dir <- "output_directory_withoutgff"  # Define your output directory
  output_dir <- create_output_dir(output_dir)
  dir.create(file.path(output_dir, "blast_result"), showWarnings = FALSE, recursive = TRUE)
  # Apply the transformations
  data_transformed_without_gff <- transform_data_without_gff(data)
  # path_data_transformed <- file.path(output_dir, paste0("data_transformed.tsv"))
  # write.table(data_transformed, file = path_data_transformed, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  # Generate the commands
  data_with_commands_without_gff <- generate_commands_without_gff(data_transformed_without_gff, output_dir)
  # Execution script
  execution_module_without_gff(data_with_commands_without_gff, output_dir)

} else if (num_cols == 5 || num_cols == 4) {
  
  source("code/moduleGff.R")

  print("with gff")
  output_dir <- "output_directory"  # Define your output directory
  output_dir <- create_output_dir(output_dir)
  dir.create(file.path(output_dir, "blast_result"), showWarnings = FALSE, recursive = TRUE)
  # Apply the transformations
  data_transformed <- transform_data(data)
  # path_data_transformed <- file.path(output_dir, paste0("data_transformed.tsv"))
  # write.table(data_transformed, file = path_data_transformed, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  # Generate the commands
  data_with_commands <- generate_commands(data_transformed, output_dir)
  # Execution script
  execution_module(data_with_commands, output_dir)

} else if (num_cols == 6) {
  
  source("code/moduleHybrid.R")

  print("Hybrid")
  output_dir <- "output_directory_hybrid"  # Define your output directory
  output_dir <- create_output_dir(output_dir)
  dir.create(file.path(output_dir, "blast_result"), showWarnings = FALSE, recursive = TRUE)
  # Apply the transformations
  data_transformed <- transform_data_hybrid(data)
  # path_data_transformed <- file.path(output_dir, paste0("data_transformed.tsv"))
  # write.table(data_transformed, file = path_data_transformed, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  # Generate the commands
  data_with_commands <- generate_commands_hybrid(data_transformed, output_dir)
  # Execution script
  execution_module_hybrid(data_with_commands, output_dir)

} else {
  stop("Unexpected number of columns in input file.")
}
