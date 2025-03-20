#!/usr/bin/env Rscript

# call all scripts
FILTER_SEQ <- "code/filter_keywords.sh"
GSCISSORS <- "code/Genomics_Scissors/gscissors.pl"
OVERLAPPINGSHAIVE <- "code/overlappingShaive.R"
BOTHBLAST <- "code/bothBlast.sh"
SPDIFFSIZE <- "code/spDiffSize.R"
GORF <- "code/gORF.sh"
ALLBLAST <- "code/allBlast.sh"
BREFINER <- "code/bRefiner.sh"
CLASSIFIER <- "code/classifier.sh"

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
  # data <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE)
  data <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

  # Determine the number of columns
  # num_cols <- ncol(data)
  num_cols <- max(colSums(data != "" & !is.na(data))) 

  # Assign names based on the number of columns
  if (num_cols == 5){
    colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length")
  } else if (num_cols == 4) {
    colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
  } else if (num_cols == 3) {
    colnames(data) <- c("fasta_file", "sequence_file", "length")
  } else if (num_cols ==2) {
    colnames(data) <- c("fasta_file", "sequence_file")
  } else {
    stop("Unexpected number of columns in input file.")
  }

  return(data)
}

# Function to transform data in feature file name with gff
transform_data <- function(data) {
  # DATA FORMATING
  data$gff_basename <- basename(data$gff_file)
  data$fasta_basename <- basename(data$fasta_file)
  data$no_gff_basename <- sub(".{4}$", "", data$gff_basename)
  data$no_fasta_basename <- sub(".{6}$", "", data$fasta_basename)
  data$keyword_sum <- paste(data$keyword1, data$keyword2, sep = "_")
  # FINISH DATA FORMATING
  
  data$filtred_name_gff <- paste0("filtered_", data$keyword_sum, "_", data$gff_basename) # filtered result of keywords
  data$non_filtred_name_gff <- paste0("non-filtered_", data$keyword2, "_", data$gff_basename) # rest of no filtered result of keywords
  data$sp_high_filtred_name_gff <- paste0("high_equal_filtered_", data$keyword_sum, "_", data$gff_basename) # filtered > or = to minimal length
  data$sp_low_filtred_name_gff <- paste0("low_filtered_", data$keyword_sum, "_", data$gff_basename) # filtered < to minimal length
  data$out_gscissors_high <- paste0("out_high_", data$keyword_sum, "_", data$no_gff_basename, ".fasta") # extract filtered >= minimal length sequence
  data$out_rest_gscissors <- paste0("out_rest_", data$keyword2, "_", data$no_gff_basename, ".fasta") # extract rest of no filtered sequence
  data$blastn_result <- paste0("blastn_", data$no_fasta_basename, ".txt") # result first BLASTN
  data$out_gscissors_high_translated <- paste0("out_high_", data$keyword_sum, "_", data$no_gff_basename, "_translated.fasta")
  # data$tblastn_result <- paste0("tblastn_", data$no_fasta_basename, ".txt") # result first PBLASTN
  data$overlappingshaive_result <- paste0("all_multigenic_family_", data$no_fasta_basename, ".tsv") # results filtered no overlap sequence
  data$overlappingshaive_result_filtered <- paste0("filtered_multigenic_family_", data$no_fasta_basename, ".tsv") # Results filtered no overlap and old sequence
  data$overlap_result_high_equal_df <- paste0("high_equal_filtered_multigenic_family_", data$no_fasta_basename, ".tsv") # Split result of no overlap > or = to minimal length
  data$overlap_result_low_df <- paste0("low_filtered_multigenic_family_", data$no_fasta_basename, ".tsv") # Split result of no overlap < to minimal length
  data$overlap_result_high_equal <- paste0("high_equal_filtered_multigenic_family_", data$no_fasta_basename, ".fasta") # extracted sequence > or = to minimal length
  data$overlap_result_low <- paste0("low_filtered_multigenic_family_", data$no_fasta_basename, ".fasta") # Extracted sequence < to minimal length
  data$gorf_result_file <- paste0("getorf_protein_", data$overlap_result_high_equal)
  data$blastp_result <- paste0("blastp_out_high_", data$keyword_sum, "_", data$no_gff_basename, ".txt")
  data$brefiner_blastp <- paste0("bRefiner_blastp_out_high_", data$keyword_sum, "_", data$no_gff_basename, ".txt")

  return(data)
}

# # Function to transform data in feature file name with any gff
# transform_data_without_gff <- function(data) {
#   data$fasta_basename <- basename(data$fasta_file)
#   data$sequence_basename <- basename(data$sequence_file)
#   data$no_fasta_basename <- sub(".{6}$", "", data$fasta_basename)
#   data$no_fasta_seq <- sub(".{6}$", "", data$sequence_basename)

#   data$blastn_result <- paste0("blastn_", data$no_fasta_basename, "_", data$no_fasta_seq, ".txt")
#   data$tblastn_result <- paste0("tblastn_", data$no_fasta_basename, "_", data$no_fasta_seq, ".txt")
#   data$overlappingshaive_result <- paste0("all_multigenic_family_", data$no_fasta_basename, "_", data$no_fasta_seq, ".tsv")
#   data$overlappingshaive_result_filtered <- paste0("filtered_multigenic_family_", data$no_fasta_basename, "_", data$no_fasta_seq, ".tsv")
#   data$overlap_result_high_equal_df <- paste0("high_equal_all_multigenic_family_", data$no_fasta_basename, "_", data$no_fasta_seq, ".tsv")
#   data$overlap_result_low_df <- paste0("low_all_multigenic_family_", data$no_fasta_basename, "_", data$no_fasta_seq, ".tsv")
#   data$overlap_result_high_equal <- paste0("high_equal_all_multigenic_family_", data$no_fasta_basename, "_", data$no_fasta_seq, ".fasta")
#   data$overlap_result_low <- paste0("low_all_multigenic_family_", data$no_fasta_basename, "_", data$no_fasta_seq, ".fasta")

#   return(data)
# }

# Function to generate commands with gff data
generate_commands <- function(data, output_dir) {
  data$filter_seq_command <- paste(FILTER_SEQ, input_file, output_dir)
  
  if ("length" %in% colnames(data)) {
    data$spdiffsize_command_first <- paste("Rscript", SPDIFFSIZE, "--gff", file.path(output_dir, data$filtred_name_gff), "--length", data$length)
  } else {
    data$spdiffsize_command_first <- paste("Rscript", SPDIFFSIZE, "--gff", file.path(output_dir, data$filtred_name_gff))
  }

  data$gscissors_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates", 
                                  file.path(output_dir, data$sp_high_filtred_name_gff), "--format gff --output", 
                                  file.path(output_dir, data$out_gscissors_high)) ##### HIGH and EQUAL #####
  # data$gscissors_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates", 
  #                                 file.path(output_dir, data$sp_low_filtred_name_gff), "--format gff --output", 
  #                                 file.path(output_dir, data$out_gscissors_high)) ##### LOW #####

  # data$gscissors_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates", 
  #                                 file.path(output_dir, data$filtred_name_gff), "--format gff --output", 
  #                                 file.path(output_dir, data$out_gscissors_high))
  data$gscissors_rest_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates", 
                                       file.path(output_dir, data$non_filtred_name_gff), "--format gff --output", 
                                       file.path(output_dir, data$out_rest_gscissors)) ###### Rest of sequence not searched ####
  
  data$bothblast_command <- paste(BOTHBLAST, file.path(output_dir, data$out_gscissors_high), file.path(output_dir, "blast_result"), data$fasta_file)
  # data$allblast_first_command <- paste(ALLBLAST, "-type", "blastn", "-qn", file.path(output_dir, data$out_gscissors_high), "-sn", data$fasta_file, "-o", file.path(output_dir, "blast_result"))
  # data$allblast_first_transeq_command <- paste(ALLBLAST, "-transeq", "-qn", file.path(output_dir, data$out_gscissors_high))
  
  data$overlappingshaive_command <- paste("Rscript", OVERLAPPINGSHAIVE, "--blast_file", file.path(output_dir, "blast_result", data$blastn_result), 
                                          "--gff_file", file.path(output_dir, data$sp_high_filtred_name_gff), "--output_dir", output_dir, "--inter", 100) ## THIS INTER OPTION NEED ESTIMATED IN PROGRAM
  
  # data$spdiffsize_command <- paste("Rscript", SPDIFFSIZE, file.path(output_dir, data$overlappingshaive_result), 8000)
  if ("length" %in% colnames(data)) {
    data$spdiffsize_command_mf <- paste("Rscript", SPDIFFSIZE, "--tsv", file.path(output_dir, data$overlappingshaive_result_filtered), "--length", data$length)
  } else {
    data$spdiffsize_command_mf <- paste("Rscript", SPDIFFSIZE, "--tsv", file.path(output_dir, data$overlappingshaive_result_filtered))
  }

  data$gscissors_high_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates",
                                  file.path(output_dir, data$overlap_result_high_equal_df), "--format", "txt", "--output",
                                  file.path(output_dir, data$overlap_result_high_equal))

  data$gscissors_low_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates",
                                       file.path(output_dir, data$overlap_result_low_df), "--format", "txt", "--output",
                                       file.path(output_dir, data$overlap_result_low))

  if ("length" %in% colnames(data)) {
    data$gorf_command <- paste("Rscript", GORF, file.path(output_dir, data$overlap_result_high_equal), output_dir, data$length)
  } else {
    data$gorf_command <- paste("Rscript", GORF, file.path(output_dir, data$overlap_result_high_equal), output_dir)
  }

  # data$allblast_blastp_command <- paste(ALLBLAST, "-type", "blastp", "-qp", file.path(output_dir, data$out_gscissors_high_translated), "-sp", data$gorf_result_file, "-o", file.path(output_dir, "blast_result"))
  
  # data$bRefiner_command <- paste(BREFINER , "-file", blastp_file, "-i", 80, "-l", data$length / 3, "-col", 1, "-unique", 1)

  data$classifier_command <- paste("Rscript", CLASSIFIER, data$overlappingshaive_result_filtered, data$brefiner_blastp)

  return(data)
}

generate_commands_without_gff <- function(data, output_dir){


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

        cat("Processing SPDIFFSIZE: ", data$spdiffsize_command[i], "\n")
        system(data$spdiffsize_command_first[i])
        # Check if SPDIFFSIZE created the expect file
        path_file_sp <- file.path(output_dir, data$sp_high_filtred_name_gff[i])
        if (!file.exists(path_file_sp)) {
          cat("Error: SPDIFFSIZE did not create the file", file.path(output_dir, data$sp_high_filtred_name_gff[i]), "and", file.path(output_dir, data$sp_low_filtred_name_gff[i]), "\n")
          next
        }   

        # change directory from where you get the data!! DEBUGGING
        cat("Processing GSCISSORS: ", data$gscissors_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gsH <- paste0(output_dir,"/", data$out_gscissors_high[i])
        system(data$gscissors_command[i])
        if (!file.exists(path_file_gsH)) {
          cat("Error: GSCISSORS did not create the file", data$out_gscissors_high[i], "\n")
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

        # cat("Processing ALLBLAST: ", data$allblast_first_command[i], "\n")
        # system(data$allblast_first_command[i])
        # # Check if ALLBLAST created the expected file
        # path_file_bb <- file.path(output_dir, "blast_result", data$blastn_result[i])
        # if (!file.exists(path_file_bb)) {
        #   cat("Error: ALLBLAST did not create the file", file.path(output_dir, "blast_result", data$blastn_result[i]), "\n")
        #   next
        # }

        # cat("Processing ALLBLAST: ", data$allblast_first_transeq_command[i], "\n")
        # system(data$allblast_first_transeq_command[i])
        # # Check if ALLBLAST created the expected file
        # path_file_first_translated <- file.path(output_dir, data$out_gscissors_high_translated[i])
        # if (!file.exists(path_file_first_translated)) {
        #   cat("Error: ALLBLAST did not create the file", file.path(output_dir, data$out_gscissors_high_translated[i]), "\n")
        #   next
        # }

        cat("Processing OVERLAPPINGSHAIVE: ", data$overlappingshaive_command[i], "\n")
        system(data$overlappingshaive_command[i])
        # Check if OVERLAPPINGSHAIVE created the expect file
        path_file_os <- file.path(output_dir, data$overlappingshaive_result[i])
        if (!file.exists(path_file_os)) {
          cat("Error: OVERLAPPINGSHAIVE did not create the file", file.path(output_dir, data$overlappingshaive_result[i]), "\n") ## debug directory
          next
        }

        cat("Processing SPDIFFSIZE: ", data$spdiffsize_command_mf[i], "\n")
        system(data$spdiffsize_command_mf[i])
        # Check if SPDIFFSIZE created the expect file
        path_file_sp_extract <- file.path(output_dir, data$overlap_result_high_equal_df[i])
        if (!file.exists(path_file_sp_extract)) {
          cat("Error: SPDIFFSIZE did not create the file", file.path(output_dir, data$overlap_result_high_equal_df[i]), "and", file.path(output_dir, data$overlap_result_low_df[i]), "\n")
          next
        }

        cat("Processing GSCISSORS: ", data$gscissors_high_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gsG <- file.path(output_dir, data$overlap_result_high_equal[i])
        system(data$gscissors_high_command[i])
        if (!file.exists(path_file_gsG)) {
          cat("Error: GSCISSORS did not create the file", data$overlap_result_high_equal[i], "\n")
          next
        }

        cat("Processing GSCISSORS: ", data$gscissors_low_command[i], "\n")
        # Check if FILTER_SEQ created the expected file
        path_file_gsP <- file.path(output_dir, data$overlap_result_low[i])
        system(data$gscissors_low_command[i])
        if (!file.exists(path_file_gsP)) {
          cat("Error: GSCISSORS did not create the file", data$overlap_result_low[i], "\n")
          next
        }

        cat("Processing GORF: ", data$gorf_command[i], "\n")
        system(data$gorf_command[i])
        # Check if SPDIFFSIZE created the expect file
        path_file_gorf <- file.path(output_dir, data$gorf_result_file[i])
        if (!file.exists(path_file_gorf)) {
          cat("Error: GORF did not create the file", file.path(output_dir, data$gorf_result_file[i]), "\n")
          next
        }

        ##################################################################################################
        # Next steps
        # think 
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

# if (num_cols == 3 || num_cols == 2) {
#   print("without gff")
# } else if (num_cols == 5 || num_cols == 4) {
#   "with gff"
# } else {
#   stop("Unexpected number of columns in input file.")
# }


# Apply the transformations
data_transformed <- transform_data(data)
# path_data_transformed <- file.path(output_dir, paste0("data_transformed.tsv"))
# write.table(data_transformed, file = path_data_transformed, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# Generate the commands
data_with_commands <- generate_commands(data_transformed, output_dir)
# Execution script
execution_module(data_with_commands, output_dir)