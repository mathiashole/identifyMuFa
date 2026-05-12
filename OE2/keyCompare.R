#!/usr/bin/env Rscript

# Charge libary
#----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Parsed arguments
#----------------------------------------------------------------------
#------------------------------------------------------------------------
parse_args_manual <- function(args) {
  # Init list to store arguments
  opts <- list(
    rib_file = NULL,
    tm_file  = NULL,
    map_file = "mapping_ids.tsv",
    output   = "coincident_domains.tsv"
  )

}
