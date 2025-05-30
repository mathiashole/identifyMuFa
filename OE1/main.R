#!/usr/bin/env Rscript

# -------------------------------------------
# EXTERNAL MODULES
# -------------------------------------------

# External scripts called during the flow
FILTER_SEQ        <- "code/filter_keywords.sh"
GSCISSORS         <- "code/Genomics_Scissors/gscissors.pl"
OVERLAPPINGSHAIVE <- "code/overlappingShaive.R"
BOTHBLAST         <- "code/bothBlast.sh"
SPDIFFSIZE        <- "code/spDiffSize.R"
GORF              <- "code/gORF.sh"
ALLBLAST          <- "code/allBlast.sh"
BREFINER          <- "code/bRefiner.sh"
CLASSIFIER        <- "code/classifier.R"
# MEANSEQ        <- "code/meanSeq.sh"

# -------------------------------------------
# AUXILIARY FUNCTIONS
# -------------------------------------------

# Create output directory if it does not exist
create_output_dir <- function(output_dir) {
  execution_path <- getwd()
  output_dir <- paste0(execution_path, "/", output_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  return(output_dir)
}

read_input <- function(input_file) {
  read.delim(input_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
}

assign_column_names <- function(data, mode) {
  n <- ncol(data)

  if (mode == "hybrid") {
    if (n == 6) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length", "sequence_file")
    } else if (n == 5) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "sequence_file")
    } else {
      stop("Unexpected format in HYBRID mode")
    }
  } else if (mode == "gff") {
    if (n == 5) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2", "length")
    } else if (n == 4) {
      colnames(data) <- c("fasta_file", "gff_file", "keyword1", "keyword2")
    } else {
      stop("Unexpected format in GFF mode")
    }
  } else if (mode == "no_gff") {
    if (n == 3) {
      colnames(data) <- c("fasta_file", "sequence_file", "length")
    } else if (n == 2) {
      colnames(data) <- c("fasta_file", "sequence_file")
    } else {
      stop("Unexpected format in NO_GFF mode")
    }
  }

  return(data)
}

# -------------------------------------------
# MANUAL ARGUMENT PARSE
# -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_arg <- NULL
mode <- NULL

# Loop through arguments
for (i in seq_along(args)) {
  if (args[i] == "--mode" || args[i] == "-m") {
    mode <- args[i + 1]
  } else if (!startsWith(args[i], "--") && is.null(input_arg)) {
    input_arg <- args[i + 1]
  }
}

# Check the number of arguments and validation
if (is.null(mode) || !(mode %in% c("hybrid", "gff", "no_gff"))) {
  stop("You must specify the mode with --mode hybrid | gff | no_gff")
}

if (is.null(input_arg)) {
  stop("You must provide the data as a file or as a tabulated line between quotes.")
}

# -------------------------------------------
# MAIN EXECUTION
# -------------------------------------------

input_file <- input_arg
data_raw <- read_input(input_arg)
data <- assign_column_names(data_raw, mode)

output_dir <- switch(mode,
  "hybrid" = "output_directory_hybrid",
  "gff" = "output_directory",
  "no_gff" = "output_directory_withoutgff"
)

source_file <- switch(mode,
  "hybrid" = "code/moduleHybrid.R",
  "gff" = "code/moduleGff.R",
  "no_gff" = "code/moduleNoGff.R"
)

source(source_file)
output_dir <- create_output_dir(output_dir)
dir.create(file.path(output_dir, "blast_result"), showWarnings = FALSE, recursive = TRUE)


if (mode == "hybrid") {
  message("Modo: HYBRID")
  transformed <- transform_data_hybrid(data)
  commands <- generate_commands_hybrid(transformed, output_dir)
  execution_module_hybrid(commands, output_dir)

} else if (mode == "gff") {
  message("Modo: GFF")
  transformed <- transform_data(data)
  commands <- generate_commands(transformed, output_dir)
  execution_module(commands, output_dir)

} else if (mode == "no_gff") {
  message("Modo: NO GFF")
  transformed <- transform_data_without_gff(data)
  commands <- generate_commands_without_gff(transformed, output_dir)
  execution_module_without_gff(commands, output_dir)
}

