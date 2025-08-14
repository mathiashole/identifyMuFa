#!/usr/bin/env Rscript

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Uso: Rscript bloutFilter.R input output EVAL_MAX SCORE_MIN COV_MIN [RATIO_MIN RATIO_MAX]")
}

infile <- args[1]
outfile <- args[2]
EVAL_MAX <- as.numeric(args[3])
SCORE_MIN <- as.numeric(args[4])
COV_MIN <- ifelse(length(args) >= 5, as.numeric(args[5]), 0)
RATIO_MIN <- ifelse(length(args) >= 6, as.numeric(args[6]), -Inf)
RATIO_MAX <- ifelse(length(args) >= 7, as.numeric(args[7]), Inf)

# init variable values
infile <- NULL
outfile <- NULL
EVAL_MAX <- NULL
SCORE_MIN <- NULL
COV_MIN <- NULL
RATIO_MIN <- -Inf
RATIO_MAX <- Inf

# Parse arguments manually
for (i in seq_along(args)) {
  if (args[i] == "--gff_file" || args[i] == "-g") {
    gff_file <- args[i + 1]
  } else if (args[i] == "--keywords" || args[i] == "-k") {
    keyword_pairs <- args[(i + 1):length(args)]
  } else if (args[i] == "--number" || args[i] == "-n") {
    number <- as.integer(args[i + 1])
  } else if (args[i] == "--strict" || args[i] == "-s") {
    strict <- TRUE
  } else if (args[i] == "--layout" || args[i] == "-l") {
    layout_id <- args[(i + 1):length(args)]
  } else if (args[i] == "--line_plot" || args[i] == "-lp") {
    line_plot <- TRUE
  } else if (args[i] == "--table" || args[i] == "-tab") {
    table_format <- args[i + 1]
  } else if (args[i] == "--colors" || args[i] == "-c") {
    colors_input <- args[(i + 1):length(args)]
  } else if (args[i] == "--palette" || args[i] == "-p") {
    palette_name <- args[i + 1]
  } else if (args[i] == "--fill_file" || args[i] == "-ff") {
    fill_file <- args[i + 1]
  } else if (args[i] == "--format" || args[i] == "-f") {
    file_format <- args[i + 1]
  } else if (args[i] == "--order_file" || args[i] == "-of") {
    order_file <- args[i + 1]
  } else if (args[i] == "--accumulated_plot" || args[i] == "-ap") {
    accumulated_plot <- TRUE
  } else if (args[i] == "--summary" || args[i] == "-sm") {
    summary_args <- TRUE
  }
}

df <- read_table2(infile, comment = "#", col_names = FALSE)

if (ncol(df) > 20) {
  # domtblout
  colnames(df)[c(3,5,13,14,16,17,18,19)] <- c(
    "tlen", "hmmlen", "Evalue_dom", "score_dom", 
    "hmmfrom", "hmmto", "alifrom", "alito"
  )
  
  df <- df %>%
    mutate(target_cov = (alito - alifrom + 1) / tlen,
           query_cov = (hmmto - hmmfrom + 1) / hmmlen,
           length_ratio = tlen / hmmlen
           ) %>%
    filter(Evalue_dom <= EVAL_MAX,
           score_dom >= SCORE_MIN,
           target_cov >= COV_MIN,
           query_cov >= COV_MIN,
           length_ratio >= RATIO_MIN,
           length_ratio <= RATIO_MAX
           )
} else {
  # tblout
  colnames(df)[c(5,6)] <- c("Evalue_full", "score_full")
  
  df <- df %>%
    filter(Evalue_full <= EVAL_MAX,
           score_full >= SCORE_MIN)
}

write_tsv(df, outfile)
