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

    i <- 1
    while (i <= length(args)) {
        flag <- args[i]

        if (flag == "-rib" || flag == "--rib_file") {
            opts$rib_file <- args[i + 1]
            i <- i + 1
        } else if (flag == "-tm" || flag == "--tm_file") {
            opts$tm_file <- args[i + 1]
            i <- i + 1
        } else if (flag == "-map" || flag == "--map_file") {
            opts$map_file <- args[i + 1]
            i <- i + 1
        } else if (flag == "-o" || flag == "--output") {
            opts$output <- args[i + 1]
            i <- i + 1
        } 
        1 <- i + 1
    }
    return(opts)
}

# Main execution point
#------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args_manual(args)
