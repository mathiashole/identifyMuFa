#!/usr/bin/env Rscript

# # Get command-line arguments
# args <- commandArgs(trailingOnly = TRUE)

# # init variable values
# blast_file <- NULL
# gff_file <- NULL

# # Parse arguments manually
# for (i in seq_along(args)) {
#   if (args[i] == "--blast_file" || args[i] == "-b") {
#     blast_file <- args[i + 1]
#   } else if (args[i] == "--gff_file" || args[i] == "-g") {
#     gff_file <- args[(i + 1):length(args)]
#     break
# args <- commandArgs(trailingOnly = TRUE)
# blast_file <- args[1]
# gff_file <- args[2]
#   }
# }

# # Validate arguments
# if (!file.exists(blast_file)) {
#   stop("The provided BLAST file does not exist.")
# }

library(readr)
library(dplyr)

# blast_data <- read.delim(blast_file, header = FALSE, sep = "\t", fill = TRUE)

blast_data <- read.delim("/home/mathias/study/maestria/scoville_data/blast_results/blastn_TriTrypDB-68_TcruziDm28c2018_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)

# Rearrange columns V9 and V10 into V15 and V16 based on condition
df <- blast_data %>%
  mutate(
    V15 = pmin(V9, V10), # Takes the smallest value between V9 and V10
    V16 = pmax(V9, V10)  # Takes the largest value between V9 and V10
  )

# Sort by V2 (group_by) and within each group by V15 and V16
df <- df %>%
  group_by(V2) %>%
  arrange(V15, V16, .by_group = TRUE) %>%
  ungroup()

# Function to eliminate overlaps and preserve the one that maximizes the difference V16 - V15
filter_non_overlapping_with_max_diff <- function(df) {
  # Split the dataframe by groups in column V2
  df_groups <- split(df, df$V2)
  
  # Initialize the final result
  result <- data.frame()
  
  # Iterate over groups
  for (group in df_groups) {
    # Sort the group by start (V15) and end (V16) positions
    group <- group[order(group$V15, group$V16), ]
    
    # Inicializar el primer grupo sin solapamientos
    non_overlapping_group <- group[1, , drop = FALSE]
    
    # Revisar las filas restantes dentro del grupo
    for (i in 2:nrow(group)) {
      last_row <- non_overlapping_group[nrow(non_overlapping_group), ]
      current_row <- group[i, ]
      
      # Si no hay solapamiento, añadir el alineamiento al grupo sin solapamientos
      if (current_row$V15 > last_row$V16) {
        non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
      } else {
        # Si hay solapamiento, comparar la diferencia entre V16 y V15
        if ((current_row$V16 - current_row$V15) > (last_row$V16 - last_row$V15)) {
          # Si el alineamiento actual es más largo, reemplazar la fila anterior
          non_overlapping_group[nrow(non_overlapping_group), ] <- current_row
        }
      }
    }
    
    # Añadir el grupo procesado al resultado final
    result <- bind_rows(result, non_overlapping_group)
  }
  
  return(result)
}

# Llamar a la función para filtrar los solapamientos
df_no_overlaps <- filter_non_overlapping_with_max_diff(df)

# Mostrar el resultado
df_no_overlaps # 1163245 1173588 vs 1163272 1174042

# Load GFF file
gff_data <- read.delim("/home/mathias/process_data/identifyMuFa/OE1/output_directory/filtered_:DGF-1_protein_coding_gene:_TriTrypDB-68_TcruziDm28c2018.gff", comment.char = "#", header = FALSE, sep = "\t")

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

filtered_data <- remove_overlaps_with_gff(df_no_overlaps, gff_data)

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
