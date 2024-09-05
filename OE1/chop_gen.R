#!/usr/bin/env Rscript

# Function to load required libraries
load_libraries <- function() {
  library(tidyverse)
}

# Function to read input file and set column names
read_input <- function(input_file) {
  data <- read_tsv(input_file, col_names = FALSE)
  colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
  return(data)
}

# Function to generate file names based on input data
generate_file_names <- function(data) {
  data %>%
    mutate(
      gff_basename = basename(gff_file),
      no_gff_basename = str_remove(gff_basename, ".{4}$"),
      keyword_sum = paste(keyword1, keyword2, sep = "_"),
      filtred_name_gff = str_c("filtered_:", keyword_sum, ":_", gff_basename),
      out_gscissors = str_c("out_:", keyword_sum, ":_", no_gff_basename, ".fasta"),
      stat_fasta_feature = str_c("stat_", keyword_sum, "_", no_gff_basename, ".tsv")
    )
}

# Function to generate commands for each step
generate_filter_seq_command <- function(input_file) {
  return(str_c("/home/usuario/BASH/chack_gff/chack_gff_v1.sh ", input_file))
}

generate_gscissors_command <- function(fasta_file, gff_file, out_gscissors) {
  return(str_c("/home/usuario/Data_Rstudio/seqExtractor/GScissors/gscissors.pl --fasta ", fasta_file, " --coordinates ", gff_file, " --format gff --output ", out_gscissors))
}

generate_fasta_feature_command <- function(out_gscissors) {
  return(str_c("/home/usuario/Data_Rstudio/chop_genome/seq_attributes/seq_attributes.R ", out_gscissors))
}

generate_distribution_command <- function(stat_fasta_feature) {
  return(str_c("/home/usuario/Data_Rstudio/statistics_of_sequence/analyze_statistics.R ", stat_fasta_feature))
}

# Function to execute a command and print the command being executed
execute_command <- function(command) {
  cat("Executing command: ", command, "\n")
  system(command)
}

# Function to process commands for each row in the data
process_commands <- function(data) {
  for (i in 1:nrow(data)) {
    cat("Processing row ", i, "\n")
    print(data[i, ])  # Print the entire row for debugging

    # Generate commands for each step    
    filter_seq_command <- generate_filter_seq_command(data$input_file[i])

    gscissors_command <- generate_gscissors_command(data$fasta_file[i], data$filtred_name_gff[i], data$out_gscissors[i])
    fasta_feature_command <- generate_fasta_feature_command(data$out_gscissors[i])
    distribution_command <- generate_distribution_command(data$stat_fasta_feature[i])

    # Execute each command and check if the expected files are created
    execute_command(filter_seq_command)

    if (!file.exists(data$filtred_name_gff[i])) {
      cat("Error: FILTER_SEQ did not create the file", data$filtred_name_gff[i], "\n")
      next
    }

    # execute_command(gscissors_command)
    # if (!file.exists(data$out_gscissors[i])) {
    #   cat("Error: GSCISSORS did not create the file", data$out_gscissors[i], "\n")
    #   next
    # }

    # execute_command(fasta_feature_command)
    # execute_command(distribution_command)

    # cat("Successfully processed", data$fasta_file[i], "and keywords:", data$keyword1[i], data$keyword2[i], "\n")
  }
}

# Main function to run the entire script
#main <- function() {
  load_libraries()

  # Check for input file argument
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("Usage: main.R <file_with_arguments.tsv>")
  }

  # Read and process the input file
  data <- read_input(args[1])
  print(data)
  data <- generate_file_names(data)
  print(data)
  process_commands(data)
#}

# Execute the main function
#main()


# #!/usr/bin/env Rscript

# # Cargar librerías necesarias
# library(tidyverse)

# # Función para procesar cada conjunto de argumentos
# process_filter_seq <- function(input_file) {
#   # Leer el archivo de entrada
#   data <- read_tsv(input_file, col_names = FALSE)
#   colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2") # Ajustar los nombres de las columnas según sea necesario

#   # Rutas completas a los scripts
#   #FILTER_SEQ <- "/home/usuario/BASH/chack_gff/chack_gff.sh"
#   FILTER_SEQ <- "/home/usuario/BASH/chack_gff/chack_gff_v1.sh"
#   GSCISSORS <- "/home/usuario/Data_Rstudio/seqExtractor/GScissors/gscissors.pl"
#   #FASTA_FEATURE <- "/home/usuario/PERL/stat_seq/stat_seq.pl"
#   SEQ_A <- "/home/usuario/Data_Rstudio/chop_genome/seq_attributes/seq_attributes.R"
#   DISTRIBUTION <- "/home/usuario/Data_Rstudio/statistics_of_sequence/analyze_statistics.R"
#   #SEQ_A <- "/home/usuario/Data_Rstudio/statistics_of_sequence/render_quarto.R"

#   # Procesar cada fila del archivo
#   data %>%
#     mutate(
#       gff_basename = basename(gff_file),
#       no_gff_basename = str_remove(gff_basename, ".{4}$"),
#       keyword_sum = paste(keyword1, keyword2, sep = "_"),
#       #filtred_name_gff = str_c("filtered_", keyword_sum, "_", gff_basename),
#       filtred_name_gff = str_c("filtered_:", keyword_sum, ":_", gff_basename),
#       # out_gscissors = str_c("out_", keyword_sum, "_", no_gff_basename, ".fasta"),
#       out_gscissors = str_c("out_:", keyword_sum, ":_", no_gff_basename, ".fasta"),
#       stat_fasta_feature = str_c("stat_", keyword_sum, "_", no_gff_basename, ".tsv")
#     ) %>%
#     rowwise() %>%
#     mutate(
#       filter_seq_command = str_c(FILTER_SEQ, " ", input_file),
#       gscissors_command = str_c(GSCISSORS, " --fasta ", fasta_file, " --coordinates ", filtred_name_gff, " --format gff --output ", out_gscissors),
#       fasta_feature_command = str_c(SEQ_A, " ", out_gscissors),
#       distribution_command = str_c(DISTRIBUTION, " ", stat_fasta_feature)
#       # fasta_feature_command = str_c(SEQ_A, " ", out_gscissors)
#     ) %>%
#     {
#       for (i in 1:nrow(.)) {
#         cat("Processing FILTER_SEQ: ", .$filter_seq_command[i], "\n")
#         system(.$filter_seq_command[i])
        
#         # Verificar si FILTER_SEQ creó el archivo esperado
#         if (!file.exists(.$filtred_name_gff[i])) {
#           cat("Error: FILTER_SEQ did not create the file", .$filtred_name_gff[i], "\n")
#           next
#         }

#         cat("Processing GSCISSORS: ", .$gscissors_command[i], "\n")
#         system(.$gscissors_command[i])
#         if (!file.exists(.$out_gscissors[i])) {
#           cat("Error: GSCISSORS did not create the file", .$out_gscissors[i], "\n")
#           next
#         }

#         cat("Processing SEQ_A: ", .$fasta_feature_command[i], "\n")
#         system(.$fasta_feature_command[i])

#         cat("Processing DISTRIBUTION: ", .$distribution_command[i], "\n")
#         system(.$distribution_command[i])

#         cat("Successfully processed", .$fasta_file[i], "and keywords:", .$keyword1[i], .$keyword2[i], "\n")
#       }
#     }
# }

# # Verificar el número de argumentos
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
#   stop("Usage: process_sequences.R <file_with_arguments.tsv>")
# }

# # Procesar archivo de entrada
# process_filter_seq(args[1])

