## FOR PARTS

library(dplyr)

blast_data <- read.delim("blastn_TriTrypDB-68_TcruziDm28c2018_Genome.txt", header = FALSE, sep = "\t", fill = TRUE)

# Reorganizar las columnas V9 y V10 en V15 y V16 según la condición
df <- blast_data %>%
  mutate(
    V15 = pmin(V9, V10), # Toma el valor menor entre V9 y V10
    V16 = pmax(V9, V10)  # Toma el valor mayor entre V9 y V10
  )

# Ordenar por V2 (group_by) y dentro de cada grupo por V15 y V16
df <- df %>%
  group_by(V2) %>%
  arrange(V15, V16, .by_group = TRUE) %>%
  ungroup()

# Elimina todas las filas duplicadas basadas en las columnas V2 y V15
df_unique <- df %>% 
  distinct(V2, V15, .keep_all = TRUE)

# Elimina todas las filas duplicadas basadas en las columnas V2 y V16
df_unique <- df_unique %>% 
  distinct(V2, V16, .keep_all = TRUE)

# Función para eliminar solapamientos y conservar el que maximiza la diferencia V16 - V15
filter_non_overlapping_with_max_diff <- function(df) {
  # Dividir el dataframe por grupos en la columna V2
  df_groups <- split(df, df$V2)
  
  # Inicializar el resultado final
  result <- data.frame()
  
  # Iterar sobre los grupos
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

### igualdad de longitud vamos por el de mayor identidad 

filter_non_overlapping_with_max_diff <- function(df) {
  # Dividir el dataframe por grupos en la columna V2
  df_groups <- split(df, df$V2)
  
  # Inicializar el resultado final
  result <- data.frame()
  
  # Iterar sobre los grupos
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
  # Dividir el dataframe por grupos en la columna V2
  df_groups <- split(df, df$V2)
  
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
