#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(dplyr)))

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file> <filter_file>")
}

# Read TSV files
input_file <- args[1]
filter_file <- args[2]

data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
filter_data <- read.table(filter_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

data <- data %>%
  mutate(V6 = paste0(V1, "_", V5))

# Modify column 5 according to matches
data$V5 <- ifelse(data$V6 %in% filter_data$V1, paste0(data$V5, "_GEN"), paste0(data$V5, "_PSEUDOGENE"))

data_out <- data[,1:5]

# Guardar el archivo modificado
output_file <- gsub("\\.tsv$", "_classified.tsv", input_file)

write.table(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Output saved to:", output_file, "\n")