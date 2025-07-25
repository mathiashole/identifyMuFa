#!/usr/bin/env bash

set -euo pipefail

# ===========================
# HELP
# ===========================
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  --aln <aln1> [aln2 ...]       Alignment files (optional if --hmm is provided)
  --hmm <hmm1> [hmm2 ...]       HMM profile files (optional if --aln is provided)
  --db <db1> [db2 ...]          Sequence database files
  --type <prot|nucl>            Type of database: 'prot' (protein) or 'nucl' (nucleotide)

Optional:
  --outdir <dir>                Output directory (default: hmm_results)
  --cpu <int>                   Number of CPUs to use (default: 1)
  --help                        Show this help message and exit
EOF
    exit 0
}

[[ $# -eq 0 || "$1" == "--help" ]] && show_help

# ===========================
# GLOBALVARIABLE
# ===========================
# Default
outdir="hmm_results"
cpu=1
# Argument parsing
aln_files=()
hmm_files=()
db_files=()
db_type=""

# ===========================
# ARGUMENT PARSING
# ===========================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --aln)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                aln_files+=("$1")
                shift
            done
            ;;
        --hmm)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                hmm_files+=("$1")
                shift
            done
            ;;
        --db)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                db_files+=("$1")
                shift
            done
            ;;
        --type)
            shift
            db_type="$1"
            shift
            ;;
        --outdir)
            shift
            outdir="$1"
            shift
            ;;
        --cpu)
            shift
            cpu="$1"
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done


# ===========================
# VALIDATIONS
# ===========================
# Validations
if [[ ${#aln_files[@]} -eq 0 && ${#hmm_files[@]} -eq 0 ]]; then
    echo "Error: Provide either --aln or --hmm files." >&2
    exit 1
fi

if [[ ${#db_files[@]} -eq 0 ]]; then
    echo "Error: No database files provided with --db" >&2
    exit 1
fi

if [[ "$db_type" != "prot" && "$db_type" != "nucl" ]]; then
    echo "Error: --type must be 'prot' or 'nucl'" >&2
    exit 1
fi

# More validations check file existence
for f in "${aln_files[@]}"; do
    [[ -f "$f" ]] || { echo "Alignment file not found: $f" >&2; exit 1; }
done

for f in "${hmm_files[@]}"; do
    [[ -f "$f" ]] || { echo "HMM file not found: $f" >&2; exit 1; }
done

for f in "${db_files[@]}"; do
    [[ -f "$f" ]] || { echo "Database file not found: $f" >&2; exit 1; }
done

# Check required tools
for tool in hmmbuild hmmsearch nhmmer; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "Required tool '$tool' not found in PATH." >&2
        exit 1
    }
done

mkdir -p "$outdir"

build_hmms() {
    local generated_hmms=()
    echo "[INFO] Building HMM profiles from alignments..."
    for aln in "${aln_files[@]}"; do
        base=$(basename "$aln")
        prefix="${base%.*}"
        hmm_out="$outdir/${prefix}.hmm"
        hmmbuild --cpu "$cpu" "$hmm_out" "$aln"
        generated_hmms+=("$hmm_out")
    done
    hmm_files+=("${generated_hmms[@]}")
}

run_search() {
    echo "[INFO] Running HMM search..."
    for hmm in "${hmm_files[@]}"; do
        hmm_base=$(basename "$hmm")
        hmm_id="${hmm_base//./_}"
        hmm_id="${hmm_id%.hmm}"

        for db in "${db_files[@]}"; do
            db_base=$(basename "$db")
            db_id="${db_base//./_}"
            db_id="${db_id%.*}"

            result_file="$outdir/${hmm_id}_vs_${db_id}.tbl"

            if [[ "$db_type" == "prot" ]]; then
                hmmsearch --cpu "$cpu" --tblout "$result_file" "$hmm" "$db"
            else
                nhmmer --cpu "$cpu" --tblout "$result_file" "$hmm" "$db"
            fi
            echo "[DONE] $result_file"
        done
    done
}

# echo "[ALL DONE] Results saved in $outdir"

[[ ${#aln_files[@]} -gt 0 ]] && build_hmms
run_search

echo "[ALL DONE] Results saved in '$outdir'"

