#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(dplyr)))

# Capturar argumentos de línea de comandos
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file> <filter_file>")
}

input_file <- "/home/mathias/study/maestria/scoville_data/filtered_multigenic_family_TriTrypDB-68_TcruziDm28c2018_Genome.tsv"
filter_file <- "/home/mathias/study/maestria/scoville_data/bRefiner_TcruziDm28c2018_Genome.tsv"

# Leer archivos TSV
input_file <- args[1]
filter_file <- args[2]

data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
filter_data <- read.table(filter_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

data <- data %>%
  mutate(V6 = paste0(V1, "_", V5))

# Modificar la columna 5 según coincidencias
data$V5 <- ifelse(data$V6 %in% filter_data$V1, paste0(data$V5, "_GEN"), paste0(data$V5, "_PSEUDOGENE"))

data_out <- data[,1:5]

# Guardar el archivo modificado
output_file <- gsub("\.tsv$", "_annotated.tsv", input_file)
write.table(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Output saved to:", output_file, "\n")