
# Function to transform data in feature file name with any gff
transform_data_without_gff <- function(data) {
  data$fasta_basename <- basename(data$fasta_file)
  data$sequence_basename <- basename(data$sequence_file)
  data$no_fasta_basename <- sub(".{6}$", "", data$fasta_basename)
  data$no_fasta_seq <- sub(".{6}$", "", data$sequence_basename)

  data$blastn_result <- paste0("blastn_", data$no_fasta_basename, ".txt")
  data$sequence_transeq <- paste0(data$no_fasta_seq, "_translated.fasta")
  data$overlappingshaive_result <- paste0("all_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlappingshaive_result_filtered <- paste0("filtered_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlap_result_high_equal_df <- paste0("high_equal_all_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlap_result_low_df <- paste0("low_all_multigenic_family_", data$no_fasta_basename, ".tsv")
  data$overlap_result_high_equal <- paste0("high_equal_all_multigenic_family_", data$no_fasta_basename, ".fasta")
  data$overlap_result_low <- paste0("low_all_multigenic_family_", data$no_fasta_basename, ".fasta")
  data$gorf_result_file <- paste0("getorf_protein_", data$overlap_result_high_equal)
  data$blastp_result <- paste0("blastp_getorf_protein_high_equal_all_multigenic_family_", data$no_fasta_basename, ".txt")
  data$brefiner_blastp <- paste0("bRefiner_", data$blastp_result)
  data$calssifier_result <- paste0("all_multigenic_family_", data$no_fasta_basename, "_classified.tsv")

  return(data)
}

# Function to generate commands with gff data
generate_commands_without_gff <- function(data, output_dir){

  # data$allblast_first_command <- paste(ALLBLAST, "-type", "blastn", "-qn", file.path(output_dir, data$sequence_file), "-sn", data$fasta_file, "-o", file.path(output_dir, "blast_result")) ## CHECK blast sequence
  data$allblast_first_command <- paste(ALLBLAST, "-type", "blastn", "-qn", data$sequence_file, "-sn", data$fasta_file, "-o", file.path(output_dir, "blast_result")) ## CHECK blast sequence
  
  # data$allblast_first_transeq_command <- paste(ALLBLAST, "-transeq", file.path(output_dir, data$sequence_file))
  data$allblast_first_transeq_command <- paste(ALLBLAST, "-transeq", data$sequence_file, "-o", output_dir)
  
  data$overlappingshaive_command <- paste("Rscript", OVERLAPPINGSHAIVE, "--blast_file", file.path(output_dir, "blast_result", data$blastn_result), 
                                          "--output_dir", output_dir, "--inter", 100) ## THIS INTER OPTION NEED ESTIMATED IN PROGRAM
  
  if ("length" %in% colnames(data)) {
    data$spdiffsize_command_mf <- paste("Rscript", SPDIFFSIZE, "--tsv", file.path(output_dir, data$overlappingshaive_result), "--length", data$length)
  } else {
    data$spdiffsize_command_mf <- paste("Rscript", SPDIFFSIZE, "--tsv", file.path(output_dir, data$overlappingshaive_result))
  }

  data$gscissors_high_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates",
                                  file.path(output_dir, data$overlap_result_high_equal_df), "--format", "txt", "--output",
                                  file.path(output_dir, data$overlap_result_high_equal))

  data$gscissors_low_command <- paste(GSCISSORS, "--fasta", data$fasta_file, "--coordinates",
                                       file.path(output_dir, data$overlap_result_low_df), "--format", "txt", "--output",
                                       file.path(output_dir, data$overlap_result_low))

  if ("length" %in% colnames(data)) {
    data$gorf_command <- paste(GORF, file.path(output_dir, data$overlap_result_high_equal), output_dir, data$length)
  } else {
    data$gorf_command <- paste(GORF, file.path(output_dir, data$overlap_result_high_equal), output_dir)
  }

  data$allblast_blastp_command <- paste(ALLBLAST, "-type", "blastp", "-qp", file.path(output_dir, data$sequence_transeq), "-sp", file.path(output_dir, data$gorf_result_file), "-o", file.path(output_dir, "blast_result"))
  
  # if ("length"%in% colnames(data)) {
  #   # data$bRefiner_command <- paste(BREFINER , "-file", file.path(output_dir, "blast_result", data$blastp_result), "-i", 80, "-l", data$length / 3, "-col", 1, "-unique", 1)
  # } else {

  # }
  data$bRefiner_command <- paste(BREFINER , "-file", file.path(output_dir, "blast_result", data$blastp_result), "-i", 80, "-l", data$length / 3, "-col", 2, "-unique", 1) ## Need mean calculated option

  data$classifier_command <- paste("Rscript", CLASSIFIER, file.path(output_dir, data$overlappingshaive_result), file.path(output_dir, "blast_result", data$brefiner_blastp))

  return(data)
}

# Function to process each set of arguments
execution_module_without_gff <- function(data, output_dir) {
      for (i in 1:nrow(data)) {
        cat("Processing ALLBLAST: ", data$allblast_first_command[i], "\n")
        system(data$allblast_first_command[i])
        # Check if ALLBLAST created the expected file
        path_file_bb <- file.path(output_dir, "blast_result", data$blastn_result[i])
        if (!file.exists(path_file_bb)) {
          cat("Error: ALLBLAST did not create the file", file.path(output_dir, "blast_result", data$blastn_result[i]), "\n")
          next
        }

        cat("Processing ALLBLAST: ", data$allblast_first_transeq_command[i], "\n")
        system(data$allblast_first_transeq_command[i])
        # Check if ALLBLAST created the expected file
        path_file_first_translated <- file.path(output_dir, data$sequence_transeq[i])
        if (!file.exists(path_file_first_translated)) {
          cat("Error: ALLBLAST did not create the file", file.path(output_dir, data$sequence_transeq[i]), "\n")
          next
        }

        cat("Processing OVERLAPPINGSHAIVE: ", data$overlappingshaive_command[i], "\n")
        system(data$overlappingshaive_command[i])
        # Check if OVERLAPPINGSHAIVE created the expect file
        path_file_os <- file.path(output_dir, data$overlappingshaive_result[i])
        if (!file.exists(path_file_os)) {
          cat("Error: OVERLAPPINGSHAIVE did not create the file", file.path(output_dir, data$overlappingshaive_result[i]), "\n")
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

        cat("Processing ALLBLAST: ", data$allblast_blastp_command[i], "\n")
        system(data$allblast_blastp_command[i])
        # Check if ALLBLAST created the expected file
        path_file_blastp <- file.path(output_dir, "blast_result", data$blastp_result[i])
        if (!file.exists(path_file_blastp)) {
          cat("Error: ALLBLAST did not create the file", file.path(output_dir, "blast_result", data$blastp_result[i]), "\n")
          next
        }

        cat("Processing BREFINER: ", data$bRefiner_command[i], "\n")
        system(data$bRefiner_command[i])
        # Check if BREFINER created the expected file
        path_file_brefiner <- file.path(output_dir, "blast_result", data$brefiner_blastp[i])
        if (!file.exists(path_file_brefiner)) {
          cat("Error: BREFINER did not create the file", file.path(output_dir, "blast_result", data$brefiner_blastp[i]), "\n")
          next
        }

        cat("Processing CLASSIFIER: ", data$classifier_command[i], "\n")
        system(data$classifier_command[i])
        # Check if CLASSIFIER created the expected file
        path_file_calssifier <- file.path(output_dir, data$calssifier_result[i])
        if (!file.exists(path_file_calssifier)) {
          cat("Error: CLASSIFIER did not create the file", data$calssifier_result[i], "\n")
          next
        }

      }
}
