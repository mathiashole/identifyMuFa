#!/bin/bash

# =============================
# Script to calculate length and GC with infoseq
# Usage: ./calculate_features.sh -type genome -out ../output/features.tsv -fasta file1.fasta file2.fasta ...
# =============================

# Global variabls
TYPE=""
OUTPUT=""
FILES=()

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -type) TYPE="$2"; shift ;;
        -out) OUTPUT="$2"; shift ;;
        -fasta) shift; while [[ "$#" -gt 0 && "$1" != -* ]]; do FILES+=("$1"); shift; done ;;
        -h|--help)
            echo "Uso: $0 -type <genome|gene> -out <output_file.tsv> -fasta <file1.fasta> [file2.fasta ...]"
            exit 0
            ;;
        *)
            echo "Unknown parametres: $1" >&2
            exit 1 ;;
    esac
    shift
done

# Validation
if [[ -z "$TYPE" || -z "$OUTPUT" || ${#FILES[@]} -eq 0 ]]; then
    echo "Error: Faltan argumentos requeridos."
    echo "Uso: $0 -type <genome|gene> -out <output_file.tsv> -fasta <file1.fasta> [file2.fasta ...]"
    exit 1
fi

# Output file
OUTPUT="features.tsv"
echo -e "file\tID\tLength\tGC_Content" > "$OUTPUT"

# Procesar cada archivo
for fasta in "${FILES[@]}"; do
    if [ ! -f "$fasta" ]; then
        echo "File not found: $fasta" >&2
        continue
    fi

    base=$(basename "$fasta")

    # Extract base name
    if [ "$TYPE" == "genome" ]; then
        extracted=$(echo "$base" | sed -E 's/.*_([^_]+)_Genome.fasta/\1/')
    elif [ "$TYPE" == "gene" ]; then
        extracted=$(echo "$base" | sed -E 's/.*_([^_]+)\.fasta/\1/')
    else
        echo "Invalid type: $TYPE (use 'genome' or 'gene')" >&2
        exit 1
    fi

    # infoseq execution
    infoseq -only -name -length -pgc "$fasta" | tail -n +2 | while read -r line; do

        ID=$(echo "$line" | awk '{print $1}')
        LENGTH=$(echo "$line" | awk '{print $2}')
        GC_RAW=$(echo "$line" | awk '{print $3}' | tr -d '%')

        # Calcular GC en decimal (dividir por 100 y redondear a 3 decimales)
        GC=$(awk "BEGIN {printf \"%.3f\", $GC_RAW/100}")

        # Parse ID
        ID=$(echo "$ID" | sed 's/=/_/g' | sed 's/[:;].*//g' | tr -d ' ')

        echo -e "${extracted}\t${ID}\t${LENGTH}\t${GC}"
    done >> "$OUTPUT"
done
