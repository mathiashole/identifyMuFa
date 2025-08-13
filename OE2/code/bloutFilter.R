#!/usr/bin/env Rscript

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Uso: Rscript filter_hmmer.R input output EVAL_MAX SCORE_MIN [COV_MIN]")
}

infile <- args[1]
outfile <- args[2]
EVAL_MAX <- as.numeric(args[3])
SCORE_MIN <- as.numeric(args[4])
COV_MIN <- ifelse(length(args) >= 5, as.numeric(args[5]), 0)
RATIO_MIN <- ifelse(length(args) >= 6, as.numeric(args[6]), -Inf)
RATIO_MAX <- ifelse(length(args) >= 7, as.numeric(args[7]), Inf)

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
