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
  --mode <all|build|search>     Choose execution mode (default: all)
  --outdir <dir>                Output directory (default: hmm_results)
  --cpu <int>                   Number of CPUs to use (default: 1)

hmmbuild options:
  --symfrac <float>             Minimum fraction of sequences with a symbol (default: none)
  --hand                        Use hand annotation (requires Stockholm format)
  --wid <float>                 Effective sequence number weighting (default: none)
  
hmmsearch output options:
  --e <float>                   E-value threshold for reporting hits
  --domE <float>                Domain E-value threshold
  --tblout <file>               Save table output
  --domtblout <file>            Save domain table output
  --pfamtblout <file>           Save pfam-style output

  --help                        Show this help message and exit
EOF
    exit 0
}

[[ $# -eq 0 || "$1" == "--help" ]] && show_help

# ===========================
# GLOBAL VARIABLE
# ===========================
# Default
outdir="hmm_results"
cpu=1
mode="all"
symfrac=""
hand=false
wid=""
wblosum=false
wpb=false
evalue=""
domevalue=""
tblout=""
domtblout=""
pfamtblout=""

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
        --aln) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do aln_files+=("$1"); shift; done ;;
        --hmm) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do hmm_files+=("$1"); shift; done ;;
        --db) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do db_files+=("$1"); shift; done ;;
        --type) shift; db_type="$1"; shift ;;
        --outdir) shift; outdir="$1"; shift ;;
        --cpu) shift; cpu="$1"; shift ;;
        --mode) shift; mode="$1"; shift ;;

        # hmmbuild options
        --symfrac) shift; symfrac="$1"; shift ;;
        --hand) hand=true; shift ;;
        --wid) shift; wid="$1"; shift ;;
        --wblosum) wblosum=true; shift ;;
        --wpb) wpb=true; shift ;;
        
        # hmmsearch options
        --e) shift; evalue="$1"; shift ;;
        --domE) shift; domevalue="$1"; shift ;;
        --tblout) shift; tblout="$1"; shift ;;
        --domtblout) shift; domtblout="$1"; shift ;;
        --pfamtblout) shift; pfamtblout="$1"; shift ;;

        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ===========================
# DIRECTORIES
# ===========================
hmms_dir="$outdir/hmms"
search_dir="$outdir/search"
logs_dir="$outdir/logs"

mkdir -p "$hmms_dir" "$search_dir" "$logs_dir"

# Central log file
logfile="$logs_dir/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$logfile") 2>&1

# Mode validation

if [[ "$mode" != "build" && "$mode" != "search" && "$mode" != "all" ]]; then
    echo "Error: --mode must be 'build', 'search', or 'all'" >&2; exit 1
fi

if [[ "$mode" != "search" && ${#aln_files[@]} -eq 0 ]]; then
    echo "Error: --aln is required in modes 'build' and 'all'" >&2; exit 1
fi

if [[ "$mode" != "build" && ${#hmm_files[@]} -eq 0 ]]; then
    echo "Error: --hmm is required in modes 'search' and 'all'" >&2; exit 1
fi

# If --wid is used but no weighting method selected, default to wblosum
if [[ -n "$wid" && "$wblosum" == false && "$wpb" == false ]]; then
    echo "[INFO] --wid detected without --wblosum or --wpb. Defaulting to --wblosum."
    wblosum=true
fi


# More validations check file existence
for f in "${aln_files[@]}"; do [[ -f "$f" ]] || { echo "Alignment file not found: $f" >&2; exit 1; }; done

for f in "${hmm_files[@]}"; do [[ -f "$f" ]] || { echo "HMM file not found: $f" >&2; exit 1; }; done

for f in "${db_files[@]}"; do [[ -f "$f" ]] || { echo "Database file not found: $f" >&2; exit 1; }; done

# Check required tools
for tool in hmmbuild hmmsearch nhmmer; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Required tool '$tool' not found in PATH." >&2; exit 1; };
done

mkdir -p "$outdir"
# mkdir -p "hmms/$outdir" "search/$outdir" "logs/$outdir"

build_hmms() {
    local generated_hmms=()
    echo "[INFO] Building HMM profiles from alignments..."
    for aln in "${aln_files[@]}"; do
        base=$(basename "$aln")
        prefix="${base%.*}"
        hmm_out="$outdir/${prefix}.hmm"

        cmd=(hmmbuild --cpu "$cpu")

        [[ "$wblosum" == true ]] && cmd+=(--wblosum)
        [[ "$wpb" == true ]] && cmd+=(--wpb)
        [[ -n "$wid" ]] && cmd+=(--wid "$wid")

        cmd+=("$hmm_out" "$aln")

        echo "[CMD] ${cmd[*]}"
        "${cmd[@]}"
        generated_hmms+=("$hmm_out")
    done
    hmm_files+=("${generated_hmms[@]}")
}

# ===========================
# HMMSEARCH
# ===========================
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

            cmd=()
            if [[ "$db_type" == "prot" ]]; then
                cmd=(hmmsearch --cpu "$cpu")
            else
                cmd=(nhmmer --cpu "$cpu")
            fi

            [[ -n "$evalue" ]] && cmd+=(--E "$evalue")
            [[ -n "$domevalue" ]] && cmd+=(--domE "$domevalue")
            [[ -n "$tblout" ]] && cmd+=(--tblout "$outdir/${hmm_id}_${db_id}.tbl")
            [[ -n "$domtblout" ]] && cmd+=(--domtblout "$outdir/${hmm_id}_${db_id}.domtbl")
            [[ -n "$pfamtblout" ]] && cmd+=(--pfamtblout "$outdir/${hmm_id}_${db_id}.pfam")

            cmd+=("$hmm" "$db")

            echo "[CMD] ${cmd[*]}"
            "${cmd[@]}"
        done
    done
}

# ===========================
# EXECUTION
# ===========================
[[ "$mode" == "build" || "$mode" == "all" ]] && build_hmms
[[ "$mode" == "search" || "$mode" == "all" ]] && run_search

echo "[ALL DONE] Results saved in '$outdir'"

