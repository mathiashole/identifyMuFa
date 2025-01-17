#!/usr/bin/env Rscript

library(dplyr)

# args <- commandArgs(trailingOnly = TRUE)
# blast_file <- args[1]

# blast_data <- read.delim(blast_file, header = FALSE, sep = "\t", fill = TRUE)

blast_data <- read.delim("blastn_TriTrypDB-68_TcruziDm28c2018_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)

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

# Remove all duplicate rows based on columns V2 and V15
df_unique <- df %>% 
  distinct(V2, V15, .keep_all = TRUE)

# Remove all duplicate rows based on columns V2 and V16
df_unique <- df_unique %>% 
  distinct(V2, V16, .keep_all = TRUE)

# Function to eliminate overlaps and preserve the one that maximizes the difference V16 - V15
filter_non_overlapping_with_max_diff <- function(df) {
  # Split the dataframe by groups in column V2
  df_groups <- split(df, df$V2)
  
  # Initialize the final result
  result <- data.frame()
  
  # Iterate over groups
  for (group in df_groups) {
    # Ordenar el grupo por las posiciones de inicio (V15) y fin (V16)
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
df_no_overlaps

# Leer el archivo GFF
gff_data <- read.delim("archivo.gff", header = FALSE, sep = "\t")

# Filtrar df_no_overlaps basado en el GFF
df_filtered <- df_no_overlaps %>%
  filter(!(
    (V1 %in% gff_data$V1) &
    ((V15 >= gff_data$V4 & V15 <= gff_data$V5) |
     (V16 >= gff_data$V4 & V16 <= gff_data$V5) |
     (V15 <= gff_data$V4 & V16 >= gff_data$V5))
  ))


### igualdad de longitud vamos por el de mayor identidad 

filter_non_overlapping_with_max_diff <- function(df) {
  # Split the dataframe by groups in column V2
  df_groups <- split(df, df$V2)
  
  # Initialize the final result
  result <- data.frame()
  
  # Iterate over groups
  for (group in df_groups) {
    # Ordenar el grupo por las posiciones de inicio (V15) y fin (V16)
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
          # Si el alineamiento actual tiene una mayor diferencia, reemplazar la fila anterior
          non_overlapping_group <- non_overlapping_group[-nrow(non_overlapping_group), ]
          non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
        } else if ((current_row$V16 - current_row$V15) == (last_row$V16 - last_row$V15)) {
          # Si la diferencia es igual, seleccionar el de mayor porcentaje de identidad (V3)
          if (current_row$V3 > last_row$V3) {
            non_overlapping_group <- non_overlapping_group[-nrow(non_overlapping_group), ]
            non_overlapping_group <- bind_rows(non_overlapping_group, current_row)
          }
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
df_no_overlaps

### extremos

filter_non_overlapping_with_extremes <- function(df) {
  # Split the dataframe by groups in column V2
  df_groups <- split(df, df$V2)
  
  # Initialize the final result
  result <- data.frame()
  
  # Iterate over groups
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
      if (current_row$V15 <= last_row$V16) {
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
