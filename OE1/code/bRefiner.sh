#!/bin/bash

# Default values
IDENTITY_MIN=0
LENGTH_MIN=0
EVALUE_MAX=1  # Default 1 (accepts all values)
COLUMNS=""
UNIQUE_COLUMNS=""

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -file) BLAST_FILE="$2"; shift ;;
        -i) IDENTITY_MIN="$2"; shift ;;
        -l) LENGTH_MIN="$2"; shift ;;
        -e) EVALUE_MAX="$2"; shift ;;
        -col) COLUMNS="$2"; shift ;;
        -unique) UNIQUE_COLUMNS="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if file is provided
if [[ -z "$BLAST_FILE" ]]; then
    echo "Error: No BLAST file provided. Use -file <filename>"
    exit 1
fi

# Filter the BLAST table using awk
awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v cols="$COLUMNS" -v unique_cols="$UNIQUE_COLUMNS" '
BEGIN { FS=OFS="\t" }
{
    if (($3 >= id_min || id_min == 0) && ($4 >= len_min || len_min == 0) && ($11 <= evalue_max || evalue_max == 1)) {
        
        # If unique columns are specified, check for duplicates
        if (unique_cols != "") {
            split(unique_cols, unique_arr, ",");
            key = "";
            for (i in unique_arr) {
                key = key $unique_arr[i] "_";
            }
            if (seen[key]++) {
                next;  # Si ya apareció, lo ignoramos
            }
        }

        # Si no se especifican columnas, imprimir toda la fila
        if (cols == "") {
            print $0;
        } else {
            split(cols, col_arr, ",");
            for (i in col_arr) {
                printf "%s\t", $col_arr[i];
            }
            printf "\n";
        }
    }
}' "$BLAST_FILE"