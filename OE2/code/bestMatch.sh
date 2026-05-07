#!/bin/bash

# --- Help / usage
usage() {
    echo "Usage: $0 -f <blast_file> [-c <colum 1|2>] [-o <output_file>]"
    echo "  -f: Input file in BLAST tabular format (outfmt 6)."
    echo "  -c: Column to use for best match (1 | 2). Default is 2."
    echo "  -o: Prefix for output file. Default is 'best_matches'."
    exit 1
}

# default values
COL=2

# --- Parse arguments
while getopts "f:c:o:" opt; do
    case $opt in
        f) INPUT=$OPTARG ;;
        c) COL=$OPTARG ;;
        o) PREFIJO=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT" ]]; then usage; fi # Check if input file is provided

# Remove extension from input file to create prefix if not provided
if [[ -z "$PREFIJO" ]]; then
    # Remove extension from input file to create prefix
    BASE_NAME=$(basename "$INPUT")
    PREFIJO="${BASE_NAME%.*}" # remove extension
    PREFIJO="${PREFIJO%.fasta}" # remove .fasta if present
    PREFIJO="${PREFIJO%.fna}"
    PREFIJO="${PREFIJO%.txt}"
fi