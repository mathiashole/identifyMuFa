#!/usr/bin/env Rscript

# UMAP Analysis with Enhanced Visualization Options
# Author: Mathias Mangino
# Version: 1.0

suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(viridis)
  library(ggrepel)
  library(patchwork)
})

# UMAP function
perform_umap_analysis <- function(data,
                                  highlight_ids = NULL,
                                  file_col = "Genome",
                                  id_col = "id",
                                  n_neighbors = 15,
                                  min_dist = 0.1,
                                  dims = 2,
                                  normalize = TRUE,
                                  point_size = 3,
                                  background_alpha = 0.1,
                                  highlight_alpha = 0.6,
                                  label_points = FALSE,
                                  plot_title = "UMAP Analysis",
                                  output_dir = "umap_results",
                                  plot_format = "png",
                                  save_plots = TRUE,
                                  col_palette = NULL) {
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Preparing data for UMAP...")
  umap_data <- data[, !(names(data) %in% c(file_col, id_col))]
  umap_data <- umap_data %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  
  clean_idx <- which(!duplicated(umap_data))
  umap_data_clean <- umap_data[clean_idx, , drop = FALSE]
  original_data_clean <- data[clean_idx, , drop = FALSE]

  if (normalize) {
    message("Normalizing data...")
    umap_data_clean <- scale(umap_data_clean)
  }

  message("Running UMAP...")
  set.seed(42)
  umap_result <- umap(umap_data_clean, n_neighbors = n_neighbors, min_dist = min_dist, n_components = dims)

  umap_scores <- as.data.frame(umap_result)
  colnames(umap_scores) <- paste0("UMAP", seq_len(dims))
  umap_scores <- cbind(original_data_clean[, c(file_col, id_col)], umap_scores)

  if (!is.null(highlight_ids)) {
    umap_scores$highlight <- ifelse(
      umap_scores[[id_col]] %in% highlight_ids,
      paste0("Highlighted_", umap_scores[[file_col]]),
      "Other"
    )
    umap_scores$alpha <- ifelse(
      grepl("Highlighted", umap_scores$highlight),
      highlight_alpha,
      background_alpha
    )
  } else {
    umap_scores$highlight <- umap_scores[[file_col]]
    umap_scores$alpha <- highlight_alpha
  }

  if (is.null(col_palette)) {
    n_colors <- length(unique(umap_scores$highlight))
    col_palette <- if (n_colors <= 10) scales::hue_pal()(n_colors) else viridis(n_colors)
    names(col_palette) <- unique(umap_scores$highlight)
    col_palette["Other"] <- "grey80"
  }

  message("Creating plot...")
  p <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = subset(umap_scores, highlight == "Other"),
               color = "grey80", alpha = background_alpha, size = point_size) +
    geom_point(aes(color = highlight, alpha = alpha), size = point_size) +
    scale_color_manual(values = col_palette) +
    scale_alpha_identity() +
    labs(
      title = plot_title,
      x = "UMAP 1",
      y = "UMAP 2",
      color = ifelse(is.null(highlight_ids), file_col, "Highlight Status")
    ) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))

  if (label_points && !is.null(highlight_ids)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(umap_scores, highlight != "Other"),
      aes(label = get(id_col)), size = 3, max.overlaps = Inf, box.padding = 0.5)
  }

#   if (save_plots) {
#     out_file <- paste0(output_prefix, ".", plot_format)
#     ggsave(out_file, p, width = 10, height = 8, dpi = 300)
#     write.csv(umap_scores, paste0(output_prefix, "_coordinates.csv"), row.names = FALSE)
#     message("Saved UMAP plot and coordinates.")
#   }

#   return(list(scores = umap_scores, plot = p))
# }
if (save_plots) {
    plot_file <- file.path(output_dir, paste0("umap_plot.", plot_format))
    coords_file <- file.path(output_dir, "umap_coordinates.csv")
    
    ggsave(plot_file, p, width = 10, height = 8, dpi = 300)
    write.csv(umap_scores, coords_file, row.names = FALSE)
    
    message("Saved UMAP plot to: ", plot_file)
    message("Saved coordinates to: ", coords_file)
  }

  return(list(scores = umap_scores, plot = p))
}

main <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Input CSV file", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = "umap_results", help = "Output prefix", metavar = "STR"),
    make_option(c("--id_col"), type = "character", default = "id", help = "ID column name", metavar = "COL"),
    make_option(c("--file_col"), type = "character", default = "Genome", help = "File/group column name", metavar = "COL"),
    make_option(c("--highlight"), type = "character", help = "File with IDs to highlight (one per line)", metavar = "FILE"),
    make_option(c("--n_neighbors"), type = "numeric", default = 15, help = "UMAP n_neighbors [default: %default]", metavar = "NUM"),
    make_option(c("--min_dist"), type = "numeric", default = 0.1, help = "UMAP min_dist [default: %default]", metavar = "NUM"),
    make_option(c("--dims"), type = "numeric", default = 2, help = "Number of UMAP dimensions", metavar = "NUM"),
    make_option(c("--point_size"), type = "numeric", default = 3, help = "Plot point size", metavar = "NUM"),
    make_option(c("--no_normalize"), action = "store_false", dest = "normalize", default = TRUE, help = "Disable normalization"),
    make_option(c("--label_points"), action = "store_true", default = FALSE, help = "Label highlighted points"),
    make_option(c("--format"), type = "character", default = "png", help = "Plot format: png, pdf, svg", metavar = "FMT")
  )

  parser <- OptionParser(option_list = option_list,
                         description = "UMAP Analysis with optional highlighting.\nExample:\n  Rscript umap_script.R --input data.csv --highlight ids.txt")
  opts <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

  if (length(commandArgs(trailingOnly = TRUE)) == 0) {
    print_help(parser)
    quit(status = 1)
  }

  if (is.null(opts$input)) {
    print_help(parser)
    stop("Input file is required.")
  }

  data <- read.csv(opts$input)
  highlight_ids <- if (!is.null(opts$highlight)) readLines(opts$highlight) else NULL

  # perform_umap_analysis(
  #   data = data,
  #   highlight_ids = highlight_ids,
  #   file_col = opts$file_col,
  #   id_col = opts$id_col,
  #   n_neighbors = opts$n_neighbors,
  #   min_dist = opts$min_dist,
  #   dims = opts$dims,
  #   normalize = opts$normalize,
  #   point_size = opts$point_size,
  #   label_points = opts$label_points,
  #   output_prefix = opts$output,
  #   plot_format = opts$format
  # )
  perform_umap_analysis(
  data = data,
  highlight_ids = highlight_ids,
  file_col = opts$file_col,
  id_col = opts$id_col,
  n_neighbors = opts$n_neighbors,
  min_dist = opts$min_dist,
  dims = opts$dims,
  normalize = opts$normalize,
  point_size = opts$point_size,
  label_points = opts$label_points,
  output_dir = opts$output,
  plot_format = opts$format
  )

}

if (!interactive()) {
  main()
}

