#!/bin/bash

# Function to display usage
usage() {
    echo "Usage:"
    echo "  - BLAST search: $0 -type <blast_type> -qp|-qn <query_fasta> -o <output_dir> -sp|-sn <genome_fasta1> [<genome_fasta2> ...]"
    echo "  - Translate only: $0 -transeq <fasta_file1> [<fasta_file2> ...]"
    echo "  -cpu <int>                 Number of threads to use (default: 1)"
    exit 1
}

# Default values
translate_only=false
blast_search=false
genomes=()
query_file=""
output_dir=""
blast_type=""
query_type=""
subject_type=""
evalue="1e-5"
threads=1

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -transeq)
            translate_only=true
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                genomes+=("$1")
                shift
            done
            ;;
        -type) blast_type="$2"; blast_search=true; shift 2;;
        -qp) query_type="protein"; query_file="$2"; shift 2;;
        -qn) query_type="nucleotide"; query_file="$2"; shift 2;;
        -sp) subject_type="protein"; genomes+=("$2"); shift 2;;
        -sn) subject_type="nucleotide"; genomes+=("$2"); shift 2;;
        -o) output_dir="$2"; shift 2;;
        -e) evalue="$2"; shift 2;;
        -cpu) threads="$2"; shift 2;;
        *) usage;;
    esac
done

echo "Query file: $query_file"
echo "Genomes: ${genomes[@]}"
echo "Output directory: $output_dir"

# Function to translate a sequence file
# translate_sequence() {
#     local input_file="$1"
#     local output_file="${input_file%.fasta}_translated.fasta"
#     transeq -clean -sequence "$input_file" -outseq "$output_file"
#     echo "$output_file"
# }
translate_sequence() {
    local input_file="$1"
    local output_dir="${2:-$(dirname "$input_file")}"
    local base_name=$(basename "$input_file")
    local output_file="${output_dir}/${base_name%.fasta}_translated.fasta"
    
    mkdir -p "$output_dir"
    transeq -clean -sequence "$input_file" -outseq "$output_file"
    echo "$output_file"
}

# If translate only mode, execute transeq
# if [[ "$translate_only" == true ]]; then
#     if [[ ${#genomes[@]} -eq 0 ]]; then
#         echo "Error: No FASTA files provided for translation."
#         exit 1
#     fi

#     for fasta in "${genomes[@]}"; do
#         if [[ ! -f "$fasta" ]]; then
#             echo "Error: File '$fasta' not found!"
#             continue
#         fi
#         translate_sequence "$fasta"
        
#     done
#     echo "Translation completed."
#     exit 0
# fi
if [[ "$translate_only" == true ]]; then
    if [[ ${#genomes[@]} -eq 0 ]]; then
        echo "Error: No FASTA files provided for translation."
        exit 1
    fi

    for fasta in "${genomes[@]}"; do
        if [[ ! -f "$fasta" ]]; then
            echo "Error: File '$fasta' not found!"
            continue
        fi
        translate_sequence "$fasta" "$output_dir"  # $output_dir
    done
    echo "Translation completed. Output saved to: ${output_dir:-$(dirname "$fasta")}"
    exit 0
fi

# If BLAST search mode, validate parameters
if [[ "$blast_search" == true ]]; then
    if [[ -z "$output_dir" || -z "$blast_type" || -z "$query_file" || -z "$query_type" || -z "$subject_type" || ${#genomes[@]} -eq 0 ]]; then
        usage
    fi

    mkdir -p "$output_dir"

    # Translate query if needed
    if [[ "$query_type" == "nucleotide" && "$blast_type" =~ ^(blastp|tblastn|blastx|psiblast)$ ]]; then
        # query_file=$(translate_sequence "$query_file")
        query_file=$(translate_sequence "$query_file" | tail -n 1)
        query_type="protein"
    fi # THIS MODULE IS VERY IMPORTANT

    # Run BLAST for each genome
    for genome in "${genomes[@]}"; do
        if [[ ! -f "$genome" ]]; then
            echo "Error: File '$genome' not found!"
            continue
        fi

        # Translate subject if needed
        if [[ "$subject_type" == "nucleotide" && "$blast_type" =~ ^(blastp|tblastn|psiblast)$ ]]; then
            genome=$(translate_sequence "$genome" | tail -n 1)
            subject_type="protein"
        fi # THIS MODULE IS VERY IMPORTANT

        # Ensure subject is in BLAST database format
        # db_prefix="${output_dir}/$(basename "$genome")"
        db_prefix="${output_dir}/$(basename "${genome%.fasta}")"

        if [[ ! -f "${db_prefix}.nhr" && ! -f "${db_prefix}.phr" ]]; then
            echo "Creating BLAST database for $genome in $output_dir..."

            if [[ "$subject_type" == "nucleotide" ]]; then
                makeblastdb -in "$genome" -dbtype nucl -out "$db_prefix"
            else
                makeblastdb -in "$genome" -dbtype prot -out "$db_prefix"
            fi
        fi

        # # output_blast="${output_dir}/$(basename "$genome")_blast_results.txt"
        # output_blast="${output_dir}/${blast_type}_$(basename "$genome").txt"
        # Define BLAST output filename (without .fasta.txt)
        output_blast="${output_dir}/${blast_type}_$(basename "${genome%.fasta}").txt"


        echo "Running BLAST ($blast_type) on $genome..."
        
        case "$blast_type" in
            blastn)
                blastn -query "$query_file" -subject "$genome" -out "$output_blast" -outfmt '6 std qlen slen' -evalue "$evalue" -num_threads "$threads"
                ;;
            blastp)
                blastp -query "$query_file" -subject "$genome" -out "$output_blast" -outfmt '6 std qlen slen' -evalue "$evalue" -num_threads "$threads"
                ;;
            tblastn)
                tblastn -query "$query_file" -subject "$genome" -out "$output_blast" -outfmt '6 std qlen slen' -evalue "$evalue" -num_threads "$threads"
                ;;
            blastx)
                blastx -query "$query_file" -subject "$genome" -out "$output_blast" -outfmt '6 std qlen slen' -evalue "$evalue" -num_threads "$threads"
                ;;
            psiblast)
                psiblast -query "$query_file" -subject "$genome" -out "$output_blast" -outfmt '6 std qlen slen' -evalue "$evalue"
                ;;
            *)
                echo "Error: Unsupported BLAST type '$blast_type'."
                exit 1
                ;;
        esac

        echo "BLAST completed for $genome. Results saved to $output_blast"
    done
    exit 0
fi

# If neither mode was triggered, show usage
usage

