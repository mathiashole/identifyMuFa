#!/usr/bin/env Rscript
# t-SNE Analysis with Enhanced Visualization Options
# Author: Mathias Mangino
# Version: 1.0

suppressPackageStartupMessages({
  library(Rtsne)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(viridis)
  library(patchwork)
})

# Enhanced t-SNE function with more options
perform_tsne_analysis <- function(data, 
                                 highlight_ids = NULL,
                                 file_col = "File", 
                                 id_col = "ID",
                                 perplexity = 30,
                                 max_iter = 1000,
                                 theta = 0.5,
                                 dims = 2,
                                 pca = TRUE,
                                 normalize = TRUE,
                                 col_palette = NULL,
                                 point_size = 3,
                                 background_alpha = 0.1,
                                 highlight_alpha = 0.6,
                                 label_points = FALSE,
                                 plot_title = "t-SNE Analysis",
                                 output_prefix = "tsne_results",
                                 save_plots = TRUE,
                                 plot_format = "png") {
  
  # Data preparation
  message("Preparing data for t-SNE...")
  tsne_data <- data[, !(names(data) %in% c(file_col, id_col))]
  
  # Convert to numeric and handle NAs
  tsne_data <- tsne_data %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  
  # Remove duplicates
  # duplicates <- duplicated(tsne_data)
  # tsne_data_clean <- tsne_data[!duplicates, ]
  # original_data_clean <- data[!duplicates, ]
  clean_idx <- which(!duplicated(tsne_data))
tsne_data_clean <- tsne_data[clean_idx, , drop = FALSE]
original_data_clean <- data[clean_idx, , drop = FALSE]

  
  # Normalize if requested
  if (normalize) {
    message("Normalizing data...")
    tsne_data_clean <- scale(tsne_data_clean)
  }
  
  # Perform t-SNE
  message("Running t-SNE...")
  set.seed(42)  # For reproducibility
  tsne_result <- Rtsne(tsne_data_clean,
                       dims = dims,
                       perplexity = perplexity,
                       theta = theta,
                       max_iter = max_iter,
                       pca = pca,
                       verbose = TRUE)
  
  # Prepare results dataframe
  tsne_scores <- as.data.frame(tsne_result$Y)
  colnames(tsne_scores) <- paste0("tSNE", 1:dims)
  
  head(tsne_scores)
  tail(tsne_scores)

  # Add metadata
  tsne_scores <- cbind(original_data_clean[, c(file_col, id_col)], tsne_scores)
  
  # Highlight logic
  if (!is.null(highlight_ids)) {
    tsne_scores$highlight <- ifelse(
      tsne_scores[[id_col]] %in% highlight_ids, 
      paste0("Highlighted_", tsne_scores[[file_col]]), 
      "Other"
    )
    
    tsne_scores$alpha <- ifelse(
      grepl("Highlighted", tsne_scores$highlight), 
      highlight_alpha, 
      background_alpha
    )
  } else {
    tsne_scores$highlight <- tsne_scores[[file_col]]
    tsne_scores$alpha <- highlight_alpha
  }
  
  # Color palette
  if (is.null(col_palette)) {
    n_colors <- length(unique(tsne_scores$highlight))
    if (n_colors <= 10) {
      col_palette <- scales::hue_pal()(n_colors)
    } else {
      col_palette <- viridis(n_colors)
    }
    names(col_palette) <- unique(tsne_scores$highlight)
    col_palette["Other"] <- "grey80"
  }
  
  # Create plot
  message("Creating visualization...")
  p <- ggplot(tsne_scores, aes(x = tSNE1, y = tSNE2)) +
    geom_point(
      data = subset(tsne_scores, highlight == "Other"),
      color = "grey80", 
      alpha = background_alpha, 
      size = point_size
    ) +
    geom_point(
      aes(color = highlight, alpha = alpha),
      size = point_size
    ) +
    scale_color_manual(values = col_palette) +
    scale_alpha_identity() +
    labs(
      title = plot_title,
      x = "t-SNE Dimension 1",
      y = "t-SNE Dimension 2",
      color = ifelse(is.null(highlight_ids), file_col, "Highlight Status")
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Add labels if requested
  if (label_points && !is.null(highlight_ids)) {
    p <- p + 
      ggrepel::geom_text_repel(
        data = subset(tsne_scores, highlight != "Other"),
        aes(label = get(id_col)),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.5
      )
  }
  
  # Save results
  if (save_plots) {
    output_file <- paste0(output_prefix, ".", plot_format)
    message("Saving plot to: ", output_file)
    ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
    
    # Save data
    write.csv(tsne_scores, paste0(output_prefix, "_coordinates.csv"), row.names = FALSE)
    message("Saved t-SNE coordinates to: ", paste0(output_prefix, "_coordinates.csv"))
  }
  
  # Return results
  return(list(
    tsne_result = tsne_result,
    tsne_scores = tsne_scores,
    tsne_plot = p,
    parameters = list(
      perplexity = perplexity,
      max_iter = max_iter,
      theta = theta,
      dims = dims,
      pca = pca,
      normalize = normalize
    )
  ))
}

# Function to run per-genome analysis
run_per_genome_analysis <- function(data, 
                                   genomes,
                                   highlight_ids,
                                   id_col = "id",
                                   output_dir = "genome_tsne",
                                   ...) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  results <- list()
  
  for (genome in genomes) {
    message("\nProcessing genome: ", genome)
    genome_data <- subset(data, Genome == genome)
    
    # Run t-SNE
    result <- perform_tsne_analysis(
      data = genome_data,
      highlight_ids = highlight_ids,
      file_col = "Genome",
      id_col = id_col,
      output_prefix = file.path(output_dir, paste0("tsne_", genome)),
      plot_title = paste("t-SNE for", genome),
      ...
    )
    
    results[[genome]] <- result
  }
  
  return(results)
}

# Command-line interface
main <- function() {
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input CSV file path", metavar="FILE"),
    make_option(c("-o", "--output"), type="character", default="tsne_results",
                help="Output directory prefix", metavar="DIR"),
    make_option(c("--id_col"), type="character", default="id",
                help="Column name for IDs", metavar="COLUMN"),
    make_option(c("--file_col"), type="character", default="Genome",
                help="Column name for file/genome grouping", metavar="COLUMN"),
    make_option(c("--highlight"), type="character", default=NULL,
                help="File with IDs to highlight (one per line)", metavar="FILE"),
    make_option(c("--perplexity"), type="numeric", default=30,
                help="t-SNE perplexity parameter", metavar="NUM"),
    make_option(c("--max_iter"), type="numeric", default=1000,
                help="Maximum t-SNE iterations", metavar="NUM"),
    make_option(c("--dims"), type="numeric", default=2,
                help="Output dimensions", metavar="NUM"),
    make_option(c("--point_size"), type="numeric", default=3,
                help="Point size in plot", metavar="NUM"),
    make_option(c("--no_pca"), action="store_false", dest="pca", default=TRUE,
                help="Disable PCA preprocessing"),
    make_option(c("--no_normalize"), action="store_false", dest="normalize", default=TRUE,
                help="Disable data normalization"),
    make_option(c("--label_points"), action="store_true", default=FALSE,
                help="Label highlighted points with their IDs"),
    make_option(c("--format"), type="character", default="png",
                help="Output plot format (png, pdf, svg)", metavar="FORMAT")
  )
  
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file must be specified", call.=FALSE)
  }
  
  # Read data
  message("Reading data from: ", opt$input)
  data <- read.csv(opt$input)
  
  # Read highlight IDs if provided
  highlight_ids <- NULL
  if (!is.null(opt$highlight)) {
    highlight_ids <- readLines(opt$highlight)
    message("Loaded ", length(highlight_ids), " IDs to highlight")
  }
  
  # Get unique genomes if doing per-genome analysis
  genomes <- unique(data[[opt$file_col]])
  
  # Run analysis
  if (length(genomes) > 1) {
    message("Detected multiple genomes (", length(genomes), "), running per-genome analysis")
    results <- run_per_genome_analysis(
      data = data,
      genomes = genomes,
      highlight_ids = highlight_ids,
      id_col = opt$id_col,
      output_dir = opt$output,
      perplexity = opt$perplexity,
      max_iter = opt$max_iter,
      dims = opt$dims,
      pca = opt$pca,
      normalize = opt$normalize,
      point_size = opt$point_size,
      label_points = opt$label_points,
      plot_format = opt$format
    )
  } else {
    message("Running single t-SNE analysis")
    results <- perform_tsne_analysis(
      data = data,
      highlight_ids = highlight_ids,
      file_col = opt$file_col,
      id_col = opt$id_col,
      perplexity = opt$perplexity,
      max_iter = opt$max_iter,
      dims = opt$dims,
      pca = opt$pca,
      normalize = opt$normalize,
      point_size = opt$point_size,
      label_points = opt$label_points,
      output_prefix = opt$output,
      plot_format = opt$format
    )
  }
  
  message("\nAnalysis complete!")
}

if (!interactive()) {
  main()
}

