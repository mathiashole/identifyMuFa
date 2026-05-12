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

# Validation of required arguments
if (is.null(opts$rib_file) || is.null(opts$tm_file)) {
    cat("Usage: ./keyCompare.R --rib_file and --tm_file are required.\n")
    quit(status = 1)
}

# 1. Charge data
#----------------------------------------------------------------------
# read the tables, asume they have tsv format
# Rename columns to have a common key for joining

clean_tm_domains <- function(df) {

    temp_df < df
    colnames(temp_df)[3] <- "st"
    colnames(temp_df)[4] <- "en"
    colnames(temp_df)[5] <- "ev"

    temp_df %>%
        distinct() %>%
        group_by(file, target_name) %>%
        # sort by start position
        arrange(st, .by_group = TRUE) %>%
        summarise(
            tlen = first(tlen),
            gap = if(n() > 1) max(lead(st) - en, na.rm = TRUE) else 0,
            # Logical priority a domain >= 370, if not, fushion
            start_tm = if(any((en - st) >= 370)) { 
                st[which((en - st) >= 370 & ev == min(ev[(en - st) >= 370]))[1]]
            } else {}
        )

}