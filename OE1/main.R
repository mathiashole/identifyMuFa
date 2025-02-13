#!/usr/bin/env Rscript

# call all scripts
FILTER_SEQ <- "code/filter_keywords.sh"
GSCISSORS <- "code/Genomics_Scissors/gscissors.pl"
OVERLAPPINGSHAIVE <- "code/overlappingShaive.R"
BOTHBLAST <- "code/bothBlast.sh"
SPDIFFSIZE <- "code/spDiffSize.R"

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
  # data <- read_tsv(input_file, col_names = FALSE)
  data <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
  return(data)
} # if you put only genome and sequence to search g and ps, this script need another execution option

# Function to transform data in feature file name with gff
transform_data <- function(data) {
  data$gff_basename <- basename(data$gff_file)
  data$fasta_basename <- basename(data$fasta_file)
  data$no_gff_basename <- sub(".{4}$", "", data$gff_basename)
  data$no_fasta_basename <- sub(".{6}$", "", data$fasta_basename)
  data$keyword_sum <- paste(data$keyword1, data$keyword2, sep = "_")
  
  data$filtred_name_gff <- paste0("filtered_", data$keyword_sum, "_", data$gff_basename)
  data$non_filtred_name_gff <- paste0("non-filtered_", data$keyword2, "_", data$gff_basename)
  data$out_gscissors <- paste0("out_", data$keyword_sum, "_", data$no_gff_basename, ".fasta")
  data$out_rest_gscissors <- paste0("out_rest", data$keyword2, "_", data$no_gff_basename, ".fasta")
  data$blastn_result <- paste0("blastn_", data$no_fasta_basename, ".txt")
  data$tblastn_result <- paste0("tblastn_", data$no_fasta_basename, ".txt")
  data$overlappingshaive_result <- paste0("all_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlappingshaive_result_filtered <- paste0("filtered_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlap_result_gene_df <- paste0("gene_all_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlap_result_pseudogene_df <- paste0("pseudogene_all_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlap_result_gene <- paste0("gene_all_multigenic_family_", data$no_fasta_basename, ".fasta")
  data$overlap_result_pseudogene <- paste0("pseudogene_all_multigenic_family_", data$no_fasta_basename, ".fasta")
  
  return(data)
}

# Function to generate commands
generate_commands <- function(data, output_dir) {
  data$filter_seq_command <- paste(FILTER_SEQ, input_file, output_dir)
  data$gscissors_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates", 
                                  file.path(output_dir, data$filtred_name_gff), "--format gff --output", 
                                  file.path(output_dir, data$out_gscissors))
  
  data$gscissors_rest_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates", 
                                       file.path(output_dir, data$non_filtred_name_gff), "--format gff --output", 
                                       file.path(output_dir, data$out_rest_gscissors))
  
  data$bothblast_command <- paste(BOTHBLAST, file.path(output_dir, data$out_gscissors), file.path(output_dir, "blast_result"), data$fasta_file)
  data$overlappingshaive_command <- paste("Rscript", OVERLAPPINGSHAIVE, "--blast_file", file.path(output_dir, "blast_result", data$blastn_result), 
                                          "--gff_file", file.path(output_dir, data$filtred_name_gff), "--output_dir", output_dir, "--inter", 100)
  
  data$spdiffsize_command <- paste("Rscript", SPDIFFSIZE, file.path(output_dir, data$overlappingshaive_result), 8000)

  data$gscissors_gene_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates",
                                  file.path(output_dir, data$overlap_result_gene_df), "--format", "txt", "--output",
                                  file.path(output_dir, data$overlap_result_gene))

  data$gscissors_pseudogene_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates",
                                       file.path(output_dir, data$overlap_result_pseudogene_df), "--format", "txt", "--output",
                                       file.path(output_dir, data$overlap_result_pseudogene))


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

        cat("Processing BOTHBLAST: ", data$bothblast_command[i], "\n")
        system(data$bothblast_command[i])
        # Check if BOTHBLAST created the expected file
        path_file_bb <- file.path(output_dir, "blast_result", data$blastn_result[i])
        if (!file.exists(path_file_bb)) {
          cat("Error: BOTHBLAST did not create the file", file.path(output_dir, "blast_result", data$blastn_result[i]), "\n")
          next
        }

        cat("Processing OVERLAPPINGSHAIVE: ", data$overlappingshaive_command[i], "\n")
        system(data$overlappingshaive_command[i])
        # Check if OVERLAPPINGSHAIVE created the expect file
        path_file_os <- file.path(output_dir, data$overlappingshaive_result[i])
        if (!file.exists(path_file_os)) {
          cat("Error: OVERLAPPINGSHAIVE did not create the file", file.path(output_dir, data$overlappingshaive_result[i]), "\n") ## debug directory
          next
        }

        cat("Processing SPDIFFSIZE: ", data$spdiffsize_command[i], "\n")
        system(data$spdiffsize_command[i])
        # Check if SPDIFFSIZE created the expect file
        path_file_sp <- file.path(output_dir, data$overlap_result_gene_df[i])
        if (!file.exists(path_file_sp)) {
          cat("Error: OVERLAPPINGSHAIVE did not create the file", file.path(output_dir, data$overlap_result_gene_df[i]), "and", file.path(output_dir, data$overlap_result_pseudogene_df[i]), "\n")
          next
        }   

        cat("Processing GSCISSORS: ", data$gscissors_gene_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gsG <- file.path(output_dir, data$overlap_result_gene[i])
        system(data$gscissors_gene_command[i])
        if (!file.exists(path_file_gsG)) {
          cat("Error: GSCISSORS did not create the file", data$overlap_result_gene[i], "\n")
          next
        }

        cat("Processing GSCISSORS: ", data$gscissors_pseudogene_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gsP <- file.path(output_dir, data$overlap_result_pseudogene[i])
        system(data$gscissors_pseudogene_command[i])
        if (!file.exists(path_file_gsP)) {
          cat("Error: GSCISSORS did not create the file", data$overlap_result_pseudogene[i], "\n")
          next
        }

        ##################################################################################################
        # Next steps
        # parsing translated sequence to obtein gene and pseudogene.
        # then think 
        ##################################################################################################

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
dir.create(file.path(output_dir, "blast_result"), showWarnings = FALSE, recursive = TRUE)

data <- read_input(input_file)

# Apply the transformations
data_transformed <- transform_data(data)
# path_data_transformed <- file.path(output_dir, paste0("data_transformed.tsv"))
# write.table(data_transformed, file = path_data_transformed, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# Generate the commands
data_with_commands <- generate_commands(data_transformed, output_dir)
# Execution script
execution_module(data_with_commands, output_dir)