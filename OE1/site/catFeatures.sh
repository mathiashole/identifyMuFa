#!/bin/bash

# =============================
# Script to calculate length and GC with infoseq
# Uso:
#   ./calculate_features.sh -type genome -out ../output/features.tsv -fasta file1.fasta file2.fasta ...
# =============================

# Default variables
TYPE=""
OUTPUT=""
FILES=()


while [[ "$#" -gt 0 ]]; do
    case $1 in
        -type) TYPE="$2"; shift ;;
        -out) OUTPUT="$2"; shift ;;
        -fasta) shift; while [[ "$#" -gt 0 && "$1" != -* ]]; do FILES+=("$1"); shift; done ;;
        -h|--help)
            echo "Usage: $0 -type <genome|gene> -out <ruta_salida.tsv> -fasta <file1.fasta> [file2.fasta ...]"
            exit 0
            ;;
        *)
            echo "Unknow parameter: $1" >&2
            exit 1 ;;
    esac
    shift
done


if [[ -z "$TYPE" || -z "$OUTPUT" || ${#FILES[@]} -eq 0 ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -type <genome|gene> -out <ruta_salida.tsv> -fasta <file1.fasta> [file2.fasta ...]"
    exit 1
fi


OUTDIR=$(dirname "$OUTPUT")
mkdir -p "$OUTDIR"


echo -e "file\tID\tLength\tGC_Content" > "$OUTPUT"


for fasta in "${FILES[@]}"; do
    if [ ! -f "$fasta" ]; then
        echo "File not found: $fasta" >&2
        continue
    fi

    base=$(basename "$fasta")


    if [ "$TYPE" == "genome" ]; then
        extracted=$(echo "$base" | sed -E 's/.*_([^_]+)_Genome.fasta/\1/')
    elif [ "$TYPE" == "gene" ]; then
        extracted=$(echo "$base" | sed -E 's/.*_([^_]+)\.fasta/\1/')
    else
        echo "invalid type: $TYPE (usa 'genome' o 'gene')" >&2
        exit 1
    fi

    
    infoseq -only -name -length -pgc "$fasta" | tail -n +2 | while read -r line; do
        ID=$(echo "$line" | awk '{print $1}' | sed 's/=/_/g' | sed 's/[:;].*//g' | tr -d ' ')
        # ID=$(echo "$line" | awk '{print $1}')
        LENGTH=$(echo "$line" | awk '{print $2}')
        GC_RAW=$(echo "$line" | awk '{print $3}' | tr -d '%')
        # echo "GC_RAW (original): $line" 
        GC=$(awk "BEGIN {printf \"%.3f\", $GC_RAW/100}")
        # GC=$(awk -v gc="$GC_RAW" 'BEGIN {printf "%.3f", gc/100}')
        echo -e "${extracted}\t${ID}\t${LENGTH}\t${GC}"
        # echo "$line"
    done >> "$OUTPUT"
done

echo "Archivo generado: $OUTPUT"
