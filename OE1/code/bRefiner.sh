#!/bin/bash

# Default values
IDENTITY_MIN=0
LENGTH_MIN=0
EVALUE_MAX=1  # Default 1 (accepts all values)
UNIQ_SORT=false  # By default, do NOT apply sort | uniq
COLUMN_TO_PRINT="ALL"  # Default: print all columns
SELF_FILTER=true  # By default, filter self-hits

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -file) BLAST_FILE="$2"; shift ;;
        -i) IDENTITY_MIN="$2"; shift ;;
        -l) LENGTH_MIN="$2"; shift ;;
        -e) EVALUE_MAX="$2"; shift ;;  # Add evalue filter
        -col) COLUMN_TO_PRINT="$2"; shift ;;
        -uniq) UNIQ_SORT=true ;;  # Apply sorting and unique filtering
        -self) SELF_FILTER="$2"; shift ;;  # Add self-hits filter option
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$BLAST_FILE" ]] || [[ -z "$LENGTH_MIN" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -file <filename> -l <length_min> [-e <evalue_max>] [-col <column_to_print>] [-i <identity_min>] [-uniq] [-self <true|false>]"
    exit 1
fi

# Get the path and base name of the file
INPUT_DIR=$(dirname "$BLAST_FILE")
INPUT_NAME=$(basename "$BLAST_FILE")

# Construct the output file name
OUTPUT_FILE="$INPUT_DIR/bRefiner_$INPUT_NAME"

# Apply filters
if [ "$UNIQ_SORT" = true ]; then
    awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v col="$COLUMN_TO_PRINT" '
    BEGIN { FS=OFS="\t" }
    ($3 >= id_min) && ($4 >= len_min) && ($11 <= evalue_max) {
        if (col == "ALL") print $0;
        else print $col;
    }' "$BLAST_FILE" | sort | uniq > "$OUTPUT_FILE"
else
    awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v col="$COLUMN_TO_PRINT" '
    BEGIN { FS=OFS="\t" }
    ($3 >= id_min) && ($4 >= len_min) && ($11 <= evalue_max) {
        if (col == "ALL") print $0;
        else print $col;
    }' "$BLAST_FILE" > "$OUTPUT_FILE"
fi


echo "Blast Refiner results saved to: $OUTPUT_FILE"