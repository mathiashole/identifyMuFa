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
