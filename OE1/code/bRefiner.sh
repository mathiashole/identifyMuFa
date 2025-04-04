#!/bin/bash

# Default values
IDENTITY_MIN=0
LENGTH_MIN=0
EVALUE_MAX=1  # Default 1 (accepts all values)
UNIQ_SORT=false  # By default, do NOT apply sort | uniq
COLUMN_TO_PRINT="ALL"  # Default: print all columns

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -file) BLAST_FILE="$2"; shift ;;
        -i) IDENTITY_MIN="$2"; shift ;;
        -l) LENGTH_MIN="$2"; shift ;;
        -e) EVALUE_MAX="$2"; shift ;;  # Add evalue filter
        -col) COLUMN_TO_PRINT="$2"; shift ;;
        -uniq) UNIQ_SORT=true ;;  # Apply sorting and unique filtering
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$BLAST_FILE" ]] || [[ -z "$LENGTH_MIN" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -file <filename> -l <length_min> [-e <evalue_max>] [-col <column_to_print>] [-i <identity_min>] [-uniq]"
    exit 1
fi

# Get the path and base name of the file
INPUT_DIR=$(dirname "$BLAST_FILE")
INPUT_NAME=$(basename "$BLAST_FILE")

# Construct the output file name
OUTPUT_FILE="$INPUT_DIR/bRefiner_$INPUT_NAME"

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
                next;  # If it already appeared, we ignore it.
            }
        }

        # If no columns are specified, print the entire row
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
}' "$BLAST_FILE" > "$OUTPUT_FILE"


echo "Blast Refiner results saved to: $OUTPUT_FILE"