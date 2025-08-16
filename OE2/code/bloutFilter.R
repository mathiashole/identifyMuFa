#!/usr/bin/env Rscript

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Uso: Rscript bloutFilter.R input output EVAL_MAX SCORE_MIN COV_MIN [RATIO_MIN RATIO_MAX]")
}

# infile <- args[1]
# outfile <- args[2]
# EVAL_MAX <- as.numeric(args[3])
# SCORE_MIN <- as.numeric(args[4])
# COV_MIN <- ifelse(length(args) >= 5, as.numeric(args[5]), 0)
# RATIO_MIN <- ifelse(length(args) >= 6, as.numeric(args[6]), -Inf)
# RATIO_MAX <- ifelse(length(args) >= 7, as.numeric(args[7]), Inf)

# init variable values
infile <- NULL
outfile <- NULL
EVAL_MAX <- NULL
SCORE_MIN <- NULL
COV_MIN <- 0
RATIO_MIN <- -Inf
RATIO_MAX <- Inf

# Parse arguments manually
for (i in seq_along(args)) {
  if (args[i] == "--input" || args[i] == "-i") {
    infile <- args[i + 1]
  } else if (args[i] == "--output" || args[i] == "-o") {
    outfile <- args[(i + 1):length(args)]
  } else if (args[i] == "--evalue" || args[i] == "-e") {
    EVAL_MAX <- as.integer(args[i + 1])
  } else if (args[i] == "--score_dom" || args[i] == "-sdom") {
    SCORE_MIN <- args[(i + 1):length(args)]
  } else if (args[i] == "--coverage_domain" || args[i] == "-cdom") {
    COV_MIN <- args[i + 1]
  } else if (args[i] == "--coverage_min" || args[i] == "-cmin") {
    RATIO_MIN <- args[(i + 1):length(args)]
  } else if (args[i] == "--coverage_max" || args[i] == "-cmax") {
    RATIO_MAX <- args[i + 1]
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
