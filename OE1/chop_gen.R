#!/usr/bin/env Rscript


# call all scripts
#FILTER_SEQ <- "/home/usuario/BASH/chack_gff/chack_gff.sh"
FILTER_SEQ <- "/home/usuario/BASH/chack_gff/chack_gff_v1.sh"
GSCISSORS <- "/home/usuario/Data_Rstudio/seqExtractor/GScissors/gscissors.pl"
#FASTA_FEATURE <- "/home/usuario/PERL/stat_seq/stat_seq.pl"
SEQ_A <- "/home/usuario/Data_Rstudio/chop_genome/seq_attributes/seq_attributes.R"
DISTRIBUTION <- "/home/usuario/Data_Rstudio/statistics_of_sequence/analyze_statistics.R"
#SEQ_A <- "/home/usuario/Data_Rstudio/statistics_of_sequence/render_quarto.R"

# Cargar librerías necesarias
library(tidyverse)

# Create output directory if it does not exist
create_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
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
    out_gscissors = str_c("out_:", keyword_sum, ":_", no_gff_basename, ".fasta"),
    stat_fasta_feature = str_c("stat_", keyword_sum, "_", no_gff_basename, ".tsv")
  )
  
  # Returns the transformed DataFrame
  return(data)
}

generate_commands <- function(data) {
  # we use rowwise to apply the function to each row individually
  data <- rowwise(data) %>%
    mutate(
      filter_seq_command = str_c(FILTER_SEQ, " ", input_file),
      gscissors_command = str_c(GSCISSORS, " --fasta ", fasta_file, " --coordinates ", filtred_name_gff, " --format gff --output ", out_gscissors),
      fasta_feature_command = str_c(SEQ_A, " ", out_gscissors),
      distribution_command = str_c(DISTRIBUTION, " ", stat_fasta_feature)
    )
  
  # Return DataFrame with generated commands
  return(data)
}


# # Función para procesar cada conjunto de argumentos
# process_filter_seq <- function(input_file) {
#   # Leer el archivo de entrada
#   # data <- read_tsv(input_file, col_names = FALSE)
#   # colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2") # Ajustar los nombres de las columnas según sea necesario

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

# Verificar el número de argumentos
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: process_sequences.R <file_with_arguments.tsv>")
}

# Execution section #

input_file <- args[1]
output_dir <- "output_directory"  # Define your output directory
data <- read_input(input_file)
# Apply the transformations
data_transformed <- transform_data(data)
# Generate the commands
data_with_commands <- generate_commands(data_transformed)
# Prit debugging
print(data_with_commands)

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

