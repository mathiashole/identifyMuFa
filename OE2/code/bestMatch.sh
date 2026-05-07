#!/bin/bash

# --- Help / usage
usage() {
    echo "Usage: $0 -f <blast_file> [-c <colum 1|2>] [-o <output_file>]"
    echo "  -f: Input file in BLAST tabular format (outfmt 6)."
    echo "  -c: Column to use for best match (1 for query, 2 for subject). Default is 1 (query)."
    echo "  -o: Prefix for output file. Default is 'best_matches'."
    exit 1
}
