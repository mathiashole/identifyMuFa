#!/bin/bash


# Function to calculate the average sequence length from a FASTA file
calculate_avg_length_fasta() {
    awk '/^>/ {if (seqlen) {sum+=seqlen; count++} seqlen=0; next} {seqlen+=length($0)} END {sum+=seqlen; count++; print int(sum/count)}' "$1"
}

# Function to calculate the average feature length from a GFF file
calculate_avg_length_gff() {
    awk '$1 !~ /^#/ {sum += ($5 - $4 + 1); count++} END {if (count > 0) print int(sum/count); else print "0"}' "$1"
}

# Check if the correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 --gff <file.gff> | --fasta <file.fasta>"
    exit 1
fi

# Process command-line options
case "$1" in
    --gff)
        calculate_avg_length_gff "$2"
        ;;
    --fasta)
        calculate_avg_length_fasta "$2"
        ;;
    *)
        echo "Invalid option. Use --gff or --fasta."
        exit 1
        ;;
esac

