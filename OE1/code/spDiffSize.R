#!/usr/bin/env Rscript

# Load required library
suppressPackageStartupMessages(library(dplyr))

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# init variable values
tsv_file <- NULL
gff_file <- NULL
threshold <- NULL

# Parse arguments manually
for (i in seq_along(args)) {
  if (args[i] == "--tsv" || args[i] == "-t") {
    tsv_file <- args[i + 1]
  } else if (args[i] == "--gff" || args[i] == "-g") {
    gff_file <- args[i + 1]
  } else if (args[i] == "--length" || args[i] == "-l") {
    threshold <- as.integer(args[i + 1])
  }
}

# Validate arguments
if (!xor(!is.null(tsv_file), !is.null(gff_file))) {
  stop("Error: You must specify either --tsv or --gff, but not both.")
}

if (!is.null(tsv_file)) {
  file_type <- "tsv"
  input_file <- tsv_file
  message("Processing TSV file: ", input_file)
  df <- read.table(input_file, header = FALSE, sep = "\t", quote = "")

  if (ncol(df) < 3 || ncol(df) > 4) stop("Error: TSV must have 3 or 4 columns.")

  df <- df %>% mutate(start = as.numeric(V2), end = as.numeric(V3))

  # Only assign 'name' if there is a 4th column
  if (ncol(df) == 4) {
    df <- df %>% mutate(name = V4)
  }

} else if (!is.null(gff_file)) {  # Read GFF file
  file_type <- "gff"
  input_file <- gff_file
  message("Processing GFF file: ", input_file)
  df <- read.delim(input_file, comment.char = "#", header = FALSE, sep = "\t") # read.table(input_file, header = FALSE, sep = "\t", quote = "", comment.char = "#")

  if (ncol(df) < 9) stop("Error: GFF must have at least 9 columns.")

  df <- df %>% mutate(start = as.numeric(V4), end = as.numeric(V5))
}

# Compute absolute difference between column 2 and 3
# df <- df %>% mutate(diff_abs = abs(V2 - V3))
df <- df %>%
  mutate(diff_abs = abs(start - end))

# Set threshold automatically if not provided
if (is.null(threshold)) {
  threshold <- mean(df$diff_abs, na.rm = TRUE)
  message("No length threshold provided. Using mean difference: ", round(threshold, 2))
}

# Split based on the threshold
df_high <- df %>% filter(diff_abs >= threshold) %>% select(-diff_abs)
df_low <- df %>% filter(diff_abs < threshold) %>% select(-diff_abs)
# df_high <- df %>% filter(diff_abs >= threshold) %>% select(V1, start, end)
# df_low <- df %>% filter(diff_abs < threshold) %>% select(V1, start, end)

# Restore column order for TSV output
if (file_type == "tsv") {
  if (ncol(df) == 4) {
    df_high <- df_high %>% select(V1, start, end, name)
    df_low <- df_low %>% select(V1, start, end, name)
  } else {
    df_high <- df_high %>% select(V1, start, end)
    df_low <- df_low %>% select(V1, start, end)
  }
}

# Get the directory and filename
output_dir <- dirname(input_file)
base_name <- tools::file_path_sans_ext(basename(input_file))

# Save the resulting dataframes
# Define output filenames based on file type
if (file_type == "tsv") {
  output_high <- file.path(output_dir, paste0("high_equal_", base_name, ".tsv"))
  output_low <- file.path(output_dir, paste0("low_", base_name, ".tsv"))
  # write.table(df_high, output_high, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.table(df_low, output_low, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
} else if (file_type == "gff") {
  output_high <- file.path(output_dir, paste0("high_equal_", base_name, ".gff"))
  output_low <- file.path(output_dir, paste0("low_", base_name, ".gff"))
  # write.table(df_high, output_high, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # write.table(df_low, output_low, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

write.table(df_high, output_high, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(df_low, output_low, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Files saved:\n", output_high, "\n", output_low, "\n")
