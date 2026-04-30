#!/bin/bash

# --- Help Function ---

usage() {
    echo "Usage: $0 -f <blast_file> [-c <1|2>] [-o <output_name_prefix>]"
    echo "  -f: input file (blast tsv)."
    echo "  -c: Column to filter by (1 for Query, 2 for Subject). Default is 2 (Subject)."
    echo "  -o: Prefix for output files. If not provided, it will be derived from the input file name."
    exit 1
}

# --- Default values ---
COL=2

# --- Processing command-line arguments ---

while getopts "f:c:o:" opt; do
    case $opt in
        f) INPUT=$OPTARG ;;
        c) COL=$OPTARG ;;
        o) PREFIJO=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT" ]]; then usage; fi

# --- Clean name for output : fasta, fna, txt ---
if [[ -z "$PREFIJO" ]]; then
    # Removed the extension to create a clean output name
    BASE_NAME=$(basename "$INPUT")
    PREFIJO="${BASE_NAME%.*}" # Remove the last extension
    PREFIJO="${PREFIJO%.fasta}" # Remove common fasta extensions
    PREFIJO="${PREFIJO%.fna}"
    PREFIJO="${PREFIJO%.txt}"
fi

# --- Logic to extract best matches ---
# 1. sort -k${COL},${COL}: gruping by the chosen ID column to ensure all entries for the same ID are together.
# 2. -k3,3nr -k4,4nr: sort by identity (column 3) and then by alignment length (column 4) in descending order, so the best match is first for each ID.
# 3. awk '!seen[$COL]++': Keep only the first occurrence of each unique ID in the chosen column, which corresponds to the best match due to the sorting.
# 4. awk final: Extract the relevant columns for coordinates. If filtering by Subject (2), we take Subject ID (2) and coordinates (9 and 10). If filtering by Query (1), we take Query ID (1) and coordinates (7 and 8).
