#!/bin/bash

# Check if getorf is installed
if ! command -v getorf &> /dev/null; then
    echo "Error: getorf is not installed or not in the PATH." >&2
    exit 1
fi

# Function to calculate the average sequence length
# calculate_avg_length() {
#     awk '/^>/ {if (seqlen) {sum+=seqlen; count++} seqlen=0; next} {seqlen+=length($0)} END {if (count>0) print int(sum/count); else print 0}' "$1"
# }
calculate_avg_length() {
    awk '/^>/ {if (seqlen) {sum+=seqlen; count++} seqlen=0; next} {seqlen+=length($0)} END {sum+=seqlen; count++; print int(sum/count)}' "$1"
}


# Check arguments
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input.fasta> <output_directory> [minsize]"
    exit 1
fi

# Assign arguments
INPUT_FASTA="$1"
OUTPUT_DIR="$2"
MINSIZE="$3"

# Validate input file
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Error: Input file '$INPUT_FASTA' does not exist." >&2
    exit 1
fi

# Validate output directory
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory '$OUTPUT_DIR' does not exist." >&2
    exit 1
fi

# If minsize is not provided, calculate the average sequence length
if [[ -z "$MINSIZE" ]]; then
    MINSIZE=$(calculate_avg_length "$INPUT_FASTA")
    if [[ "$MINSIZE" -le 0 ]]; then
        echo "Error: Unable to calculate a valid minsize." >&2
        exit 1
    fi
    echo "Using calculated minsize: $MINSIZE"
fi

adjusted_minsize=$(( minsize * 90 / 100 ))

# Get the input filename without path
BASENAME=$(basename "$INPUT_FASTA")

# Construct output file name
# OUTPUT_FILE="$OUTPUT_DIR/getorf_${MINSIZE}_${BASENAME}"
# OUTPUT_FILE="$OUTPUT_DIR/getorf_${adjusted_minsize}_${BASENAME}"

# Run getorf
# getorf -sequence "$INPUT_FASTA" -outseq "$OUTPUT_FILE" -minsize "$MINSIZE"
getorf -sequence "$INPUT_FASTA" -outseq "$OUTPUT_FILE" -minsize "$adjusted_minsize"

# Check if getorf succeeded
if [[ $? -eq 0 ]]; then
    echo "getorf completed successfully. Output saved to: $OUTPUT_FILE"
else
    echo "Error: getorf failed." >&2
    exit 1
fi
