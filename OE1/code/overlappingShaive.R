#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# init variable values
blast_file <- NULL
gff_file <- NULL
table_format <- NULL
inter_seq <- 0

# Parse arguments manually
for (i in seq_along(args)) {
  if (args[i] == "--blast_file" || args[i] == "-b") {
    blast_file <- args[i + 1]
  } else if (args[i] == "--gff_file" || args[i] == "-g") {
    gff_file <- args[i + 1]
  } else if (args[i] == "--inter" || args[i] == "-i") {
    inter_seq <- as.integer(args[i + 1])
  } else if (args[i] == "--table" || args[i] == "-tab") {
    table_format <- args[i + 1]
  } else if (args[i] == "--output_dir" || args[i] == "-o") {
    output_dir <- args[i + 1]  # Get the output directory
  }
}

# Validate arguments
if (!file.exists(blast_file)) {
  stop("The provided BLAST file does not exist.")
}

if (is.null(table_format)) {
  table_format <- "tsv"
}

library(readr)
library(dplyr)

blast_data <- read.delim(blast_file, header = FALSE, sep = "\t", fill = TRUE)

# inter_seq <- 100
# blast_data <- read.delim("/home/mathias/study/maestria/scoville_data/blast_results/blastn_TriTrypDB-68_TcruziDm28c2018_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)
# blast_data <- read.delim("/home/mathias/study/maestria/scoville_data/blast_results/blastn_TriTrypDB-68_TcruziBrazilA4_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)
# blast_data <- read.delim("/home/mathias/study/maestria/scoville_data/blast_results/blastn_TriTrypDB-68_TcruziYC6_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)
# blast_data <- read.delim("/home/mathias/study/maestria/scoville_data/blast_results/blastn_TcDm25_TcruziTcDm25H1_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)


# Rearrange columns V9 and V10 into V15 and V16 based on condition
df <- blast_data %>%
  mutate(
    V15 = pmin(V9, V10), # Takes the smallest value between V9 and V10
    V16 = pmax(V9, V10),  # Takes the largest value between V9 and V10
    V17 = ifelse(V9 <= V10, "+", "-")  # "+" if V9 is the minimum, "-" if V9 is the maximum
  )

# Sort by V2 (group_by) and within each group by V15 and V16
df <- df %>%
  # group_by(V2) %>%
  group_by(V2, V17) %>% # Group by V2 and then by V17
  arrange(V15, V16, .by_group = TRUE) %>% # Sort within each subgroup
  ungroup()

filter_non_overlapping_with_extremes <- function(df) {
  # # Dividir el dataframe por grupos en la columna V2
  # df_groups <- split(df, df$V2)
  # Dividir el dataframe por grupos en las columnas V2 y V17
  df_groups <- split(df, interaction(df$V2, df$V17, drop = TRUE)) # df_no_overlaps[,c(2,9,10,15,16, 17)]
  
  # Inicializar el resultado final
  result <- data.frame()
  
  # Iterar sobre los grupos
  for (group in df_groups) {
    # Ordenar el grupo por las posiciones de inicio (V15)
    group <- group[order(group$V15), ]
    
    # Inicializar el primer grupo sin solapamientos
    non_overlapping_group <- group[1, , drop = FALSE]
    
    # Revisar las filas restantes dentro del grupo
    for (i in 2:nrow(group)) {
      last_row <- non_overlapping_group[nrow(non_overlapping_group), ]
      current_row <- group[i, ]
      
      # Si hay solapamiento, actualizar los extremos
      if (current_row$V15 <= last_row$V16 + inter_seq) {
        last_row$V15 <- min(last_row$V15, current_row$V15)
        last_row$V16 <- max(last_row$V16, current_row$V16)
        non_overlapping_group[nrow(non_overlapping_group), ] <- last_row
      } else {
        # Si no hay solapamiento, añadir el alineamiento al grupo sin solapamientos
        non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
      }
    }
    
    # Añadir el grupo procesado al resultado final
    result <- bind_rows(result, non_overlapping_group)
  }
  
  return(result)
}

# Llamar a la función para filtrar los solapamientos
df_no_overlaps <- filter_non_overlapping_with_extremes(df)

# Mostrar el resultado
df_no_overlaps

# Load GFF file
# gff_data <- read.delim("/home/mathias/process_data/identifyMuFa/OE1/output_directory/filtered_:DGF-1_protein_coding_gene:_TriTrypDB-68_TcruziDm28c2018.gff", comment.char = "#", header = FALSE, sep = "\t")
# gff_data <- read.delim("/home/mathias/process_data/identifyMuFa/OE1/output_directory/filtered_:DGF-1_protein_coding_gene:_TriTrypDB-68_TcruziBrazilA4.gff", comment.char = "#", header = FALSE, sep = "\t")
# gff_data <- read.delim("/home/mathias/process_data/identifyMuFa/OE1/output_directory/filtered_:DGF-1_protein_coding_gene:_TriTrypDB-68_TcruziYC6.gff", comment.char = "#", header = FALSE, sep = "\t")
# gff_data <- read.delim("/home/mathias/process_data/identifyMuFa/OE1/output_directory/filtered_:DGF-1_polypeptide:_TcDm25_TcruziTcDm25H1.gff", comment.char = "#", header = FALSE, sep = "\t")

remove_overlaps_with_gff <- function(dataframe, gff) {
  # Filtrar las filas del dataframe comparando con los datos del GFF
  filtered_data <- dataframe[!apply(dataframe, 1, function(row) {
    # Buscar en el GFF donde V1 coincide con V2 del dataframe
    matches <- gff[gff$V1 == row["V2"], ]
    
    # Verificar si hay solapamiento con los intervalos
    any(apply(matches, 1, function(match_row) {
      # Comparar si hay solapamiento
      overlap <- (as.numeric(row["V15"]) <= as.numeric(match_row["V5"]) &&
                  as.numeric(row["V16"]) >= as.numeric(match_row["V4"]))
      return(overlap)
    }))
  }), ]
  
  return(filtered_data)
}

if (!is.null(gff_file)) {
  filtered_data <- remove_overlaps_with_gff(df_no_overlaps, gff_data)
} else {
  filtered_data <- df_no_overlaps  # If there is no GFF, keep the data unchanged
}

# We extract the base name of the blast_file file without the extension
blast_base <- tools::file_path_sans_ext(basename(blast_file))
# Remove the "blastn_" prefix if present
blast_base <- sub("^blastn_", "", blast_base)

if (!is.null(table_format)) {
  if (table_format == "csv") {
    write_csv(df_no_overlaps, file.path(output_dir, paste0("all_multigenic_family_", blast_base, ".csv")), col_names = FALSE)
    write_csv(filtered_data, file.path(output_dir, paste0("filtered_multigenic_family_", blast_base, ".csv")), col_names = FALSE)
  } else if (table_format == "tsv") {
    write_tsv(df_no_overlaps, file.path(output_dir, paste0("all_multigenic_family_", blast_base, ".tsv")), col_names = FALSE)
    write_tsv(filtered_data, file.path(output_dir, paste0("filtered_multigenic_family_", blast_base, ".tsv")), col_names = FALSE)
  } else {
    stop("Unsupported table format specified.")
  }
  cat("Tables saved in", table_format, "format without column names.\n")
}

# ### igualdad de longitud vamos por el de mayor identidad 

# filter_non_overlapping_with_max_diff <- function(df) {
#   # Split the dataframe by groups in column V2
#   df_groups <- split(df, df$V2)
  
#   # Initialize the final result
#   result <- data.frame()
  
#   # Iterate over groups
#   for (group in df_groups) {
#     # Ordenar el grupo por las posiciones de inicio (V15) y fin (V16)
#     group <- group[order(group$V15, group$V16), ]
    
#     # Inicializar el primer grupo sin solapamientos
#     non_overlapping_group <- group[1, , drop = FALSE]
    
#     # Revisar las filas restantes dentro del grupo
#     for (i in 2:nrow(group)) {
#       last_row <- non_overlapping_group[nrow(non_overlapping_group), ]
#       current_row <- group[i, ]
      
#       # Si no hay solapamiento, añadir el alineamiento al grupo sin solapamientos
#       if (current_row$V15 > last_row$V16) {
#         non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
#       } else {
#         # Si hay solapamiento, comparar la diferencia entre V16 y V15
#         if ((current_row$V16 - current_row$V15) > (last_row$V16 - last_row$V15)) {
#           # Si el alineamiento actual tiene una mayor diferencia, reemplazar la fila anterior
#           non_overlapping_group <- non_overlapping_group[-nrow(non_overlapping_group), ]
#           non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
#         } else if ((current_row$V16 - current_row$V15) == (last_row$V16 - last_row$V15)) {
#           # Si la diferencia es igual, seleccionar el de mayor porcentaje de identidad (V3)
#           if (current_row$V3 > last_row$V3) {
#             non_overlapping_group <- non_overlapping_group[-nrow(non_overlapping_group), ]
#             non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
#           }
#         }
#       }
#     }
    
#     # Añadir el grupo procesado al resultado final
#     result <- bind_rows(result, non_overlapping_group)
#   }
  
#   return(result)
# }

# # Llamar a la función para filtrar los solapamientos
# df_no_overlaps <- filter_non_overlapping_with_max_diff(df)

# # Mostrar el resultado
# df_no_overlaps

## origianl best dif

# # Function to eliminate overlaps and preserve the one that maximizes the difference V16 - V15
# filter_non_overlapping_with_max_diff <- function(df) {
#   # Split the dataframe by groups in column V2
#   df_groups <- split(df, df$V2)
  
#   # Initialize the final result
#   result <- data.frame()
  
#   # Iterate over groups
#   for (group in df_groups) {
#     # Sort the group by start (V15) and end (V16) positions
#     group <- group[order(group$V15, group$V16), ]
    
#     # Initialize the first group without overlaps
#     non_overlapping_group <- group[1, , drop = FALSE]
    
#     # Check the remaining rows within the group
#     for (i in 2:nrow(group)) {
#       last_row <- non_overlapping_group[nrow(non_overlapping_group), ]
#       current_row <- group[i, ]
      
#       # If there is no overlap, add the alignment to the group without overlaps
#       if (current_row$V15 > last_row$V16) {
#         non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
#       } else {
#         # If there is overlap, compare the difference between V16 and V15
#         if ((current_row$V16 - current_row$V15) > (last_row$V16 - last_row$V15)) {
#           # If the current alignment is longer, replace the previous row
#           non_overlapping_group[nrow(non_overlapping_group), ] <- current_row
#         }
#       }
#     }
    
#     # Add the processed group to the final result
#     result <- bind_rows(result, non_overlapping_group)
#   }
  
#   return(result)
# }

# # Call the function to filter the overlaps
# df_no_overlaps <- filter_non_overlapping_with_max_diff(df)

# # Show the result
# df_no_overlaps # 1163245 1173588 vs 1163272 1174042