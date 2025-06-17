#!/bin/bash

# Check if getorf is installed
if ! command -v getorf &> /dev/null; then
    echo "Error: getorf is not installed or not in the PATH." >&2
    exit 1
fi

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

adjusted_minsize=$(( MINSIZE * 5 / 100 ))
# AMINOACID_MINSIZE=$(( MINSIZE / 3 ))

# Get the input filename without path
BASENAME=$(basename "$INPUT_FASTA")

# Construct output file name
# OUTPUT_FILE="$OUTPUT_DIR/getorf_${MINSIZE}_${BASENAME}"
# OUTPUT_FILE="$OUTPUT_DIR/getorf_${adjusted_minsize}_${BASENAME}"
OUTPUT_FILE="$OUTPUT_DIR/getorf_${BASENAME}"
FILTERED_OUTPUT_FILE="$OUTPUT_DIR/getorf_filtered_${BASENAME}"
OUTPUT_FILE_TRANSEQ="$OUTPUT_DIR/getorf_protein_${BASENAME}"
OUTPUT_FILTERED_FILE_TRANSEQ="$OUTPUT_DIR/getorf_filtered_protein_${BASENAME}"

# Run getorf nucleotide
getorf -sequence "$INPUT_FASTA" -outseq "$OUTPUT_FILE" -minsize "$MINSIZE" -find 3
sed -i -E 's/_1(\s*\[)/\1/' "$OUTPUT_FILE"

getorf -sequence "$INPUT_FASTA" -outseq "$FILTERED_OUTPUT_FILE" -minsize "$adjusted_minsize" -maxsize "$MINSIZE" -find 3
sed -i -E 's/_1(\s*\[)/\1/' "$FILTERED_OUTPUT_FILE"

# Run getorf aminoacid
getorf -sequence "$INPUT_FASTA" -outseq "$OUTPUT_FILE_TRANSEQ" -minsize "$MINSIZE" -find 1
sed -i -E 's/_1(\s*\[)/\1/' "$OUTPUT_FILE_TRANSEQ"

getorf -sequence "$INPUT_FASTA" -outseq "$OUTPUT_FILTERED_FILE_TRANSEQ" -minsize "$adjusted_minsize" -maxsize "$MINSIZE" -find 1
sed -i -E 's/_1(\s*\[)/\1/' "$OUTPUT_FILTERED_FILE_TRANSEQ"

# Check if getorf succeeded
if [[ $? -eq 0 ]]; then
    echo "getorf completed successfully. Output saved to: $OUTPUT_FILE"
else
    echo "Error: getorf failed." >&2
    exit 1
fi
