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
            # Calculate gap between consecutive domains
            gap = if(n() > 1) max(lead(st) - en, na.rm = TRUE) else 0,
            # Logical priority a domain >= 370, if not, fushion
            start_tm = if(any((en - st) >= 370)) { 
                st[which((en - st) >= 370 & ev == min(ev[(en - st) >= 370]))[1]]
            } else if (gap < 100) {
                min(st)
            } else {
                st[which(ev == min(ev))[1]]
            },

            end_tm = if(any((en - st) >= 370)) { 
                en[which((en - st) >= 370 & ev == min(ev[(en - st) >= 370]))[1]]
            } else if (gap < 100) {
                max(en)
            } else {
                en[which(ev == min(ev))[1]]
            },

            evalue_tm = min(ev),
            .groups = "drop"
        ) %>%
    # final result take correct name and select columns
    select(file, target_name, start_tm, end_tm, evalue_tm, tlen)
}

rib_df <- read_tsv(opts$rib_file, col_types = cols()) %>%
    rename(start_rib = start, end_rib = end, evalue_rib = e_value, tlen = tlen)

rib_df <- filter(rib_df, (end_rib - start_rib) >= mean(tlen)*0.4)

tm_df <- read_tsv(opts$tm_file, col_types = cols()) %>%
    rename(start_tm = start, end_tm = end, evalue_tm = e_value, tlen = tlen)

tm_df <- clean_tm_domains(tm_df)

tm_df <- filter(tm_df, (end_tm - start_tm) >= (tlen*0.35))

map_df <- read_tsv(opts$map_file, col_names = FALSE, col_types = cols()) %>%
    rename(file = X2, target_name = X1)

map_df <- map_df %>%
    distinct(file)

# 2. Join tables
#----------------------------------------------------------------------

# inner join rib_df and tm_df on file and target_name conserving only the rows with the same file and target_name in both tables
merged_df <- inner_join(rib_df, tm_df, by = c("file", "target_name"))

# Frequency only have both domains
freq_merged <- merged_df %>%
    group_by(file, target_name) %>%
    summarise(count = n(), .groups = "drop")

# IDs with only rib domains
only_rib_df <- anti_join(rib_df, tm_df, by = c("file", "target_name"))

# Frequency only have rib domains
freq_only_rib <- only_rib_df %>%
    group_by(file, target_name) %>%
    summarise(count = n(), .groups = "drop")

# IDs with only tm domains
only_tm_df <- anti_join(tm_df, rib_df, by = c("file", "target_name"))

# Frequency only have tm domains
freq_only_tm <- only_tm_df %>%
    group_by(file, target_name) %>%
    summarise(count = n(), .groups = "drop")

# 3. Count frequencies per specie (column 'file')
#----------------------------------------------------------------------

summary_df <- merged_df %>%
    group_by(file) %>%
    summarise(both_domains = n(),
        uniq_domain = n_distinct(target_name)
    ) %>%
    arrange(desc(both_domains))

# 4. Count domain per category
#----------------------------------------------------------------------

count_rib <- only_rib_df %>%
    group_by(file) %>%
    summarise(rib = n()) # count total rib column group by file

count_tm <- only_tm_df %>%
    group_by(file) %>%
    summarise(tm = n_distinct(target_name)) # count distinct tm in target_name and column group by file

count_both <- merged_df %>%
    group_by(file) %>%
    summarise(both = n_distinct(target_name)) # count distinct target names with both domains and column group by file

both_7rib <- merged_df %>%
    group_by(file, target_name) %>%
    summarise(n_rib = n_distinct(paste(start_rib, end_rib)), .groups = "drop") %>% # count distinct rib domains because obtein start and end coordinates
    filter(n_rib == 7) %>% # filtered of number of rib
    group_by(file) %>%
    summarise(DGF1 = n_distinct(target_name), .groups = "drop") # finaly count distinct target names group by file

# 5. Read a mapping dataframe obtein all spacies
#----------------------------------------------------------------------

summary_table <- mapp_df %>%
    full_join(count_rib, by = "file") %>% # join with the dataframe containing IDs with rib domains
    full_join(count_tm, by = "file") %>% # join with the dataframe containing IDs with tm domains
    full_join(count_both, by = "file") %>% # join with the dataframe containing IDs with both domains
    full_join(both_7rib, by = "file") %>% # join with the dataframe containing IDs with 7 rib domains
    mutate(rib = coalesce(rib, 0), tm = coalesce(tm, 0), both = coalesce(both, 0), DGF1 = coalesce(DGF1, 0)) # fill NAs with 0

sevenRibId <- merged_df %>%
    group_by(file, target_name) %>%
    summarise(n_rib = n_distinct(paste(start_rib, end_rib)), .groups = "drop") %>% # count distinct rib domains because obtein start and end coordinates
    group_by(file) %>%
    filter(n_rib == 7) %>% # filtered of number of rib