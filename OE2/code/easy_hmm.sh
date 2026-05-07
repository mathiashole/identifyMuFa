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

Modo Recomendación:
  --dry-run                     Analyze and print commands without executing (for testing)

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
tblout=false
domtblout=false
pfamtblout=false
dry_run=false

# Argument parsing
aln_files=()
hmm_files=()
db_files=()
db_type=""

# Timestamp for log files (hour_daymonthyear)
timestamp=$(date +"%H%M_%d%m%Y")

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

        --dry-run) dry_run=true; shift ;;

        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# VALIDATIONS
mkdir -p "$outdir"

# Check required tools
for tool in hmmbuild hmmsearch nhmmer infoalign bc; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Error: '$tool' no instalado." >&2; exit 1; }
done

# ALIGNMENT ANALYSIS FUNCTION

analyze_alignment() {
    local aln=$1
    local output
    
    # Trying to run infoalign, if it fails, print unknown metrics and return
    if ! output=$(infoalign -sequence "$aln" -filter -auto 2>/dev/null); then
        echo "METRICS|ID:unknown|GAPS:unknown" >&2
        return
    fi

    local metrics # metrics will hold the average identity and gap percentage
    metrics=$(echo "$output" | grep -v "^#" | grep -v "USA" | awk '
        BEGIN {sum_id=0; sum_gl=0; sum_al=0; count=0} 
        {
            # Basado en tu output:
            # $4 = AlignLen
            # $6 = GapLen
            # $7 = Ident
            if($4 > 0) {
                sum_id += $7; 
                sum_gl += $6; 
                sum_al += $4; 
                count++;
            }
        } 
        END {
            if(count > 0) {
                # Identidad media (%) y Gaps medios (%)
                printf "%.2f %.2f", (sum_id/sum_al)*100, (sum_gl/sum_al)*100
            } else {
                printf "0 0"
            }
        }')

    local avg_id=$(echo "$metrics" | cut -d' ' -f1)
    local gap_pct=$(echo "$metrics" | cut -d' ' -f2)

    echo "METRICS|ID:${avg_id}%|GAPS:${gap_pct}%" >&2

    local recs=()
    # Conditional recommendations based on metrics
    if (( $(echo "$avg_id < 25" | bc -l) )); then
        recs+=("--wblosum") # If identity is very low, use BLOSUM weighting to help with distant relationships
    elif (( $(echo "$avg_id > 85" | bc -l) )); then
        recs+=("--wpb") # If identity is very high, use position-based weighting to avoid over-representation of nearly identical sequences
    fi

    if (( $(echo "$gap_pct > 30" | bc -l) )); then
        recs+=("--symfrac 0.2") # If there are many gaps, require a higher fraction of sequences to have a symbol to include that column
    elif (( $(echo "$gap_pct < 10" | bc -l) )); then
        recs+=("--symfrac 0.6") # If there are few gaps, can be more lenient with including columns
    fi
    echo "${recs[@]}"
}

# HMM BUILDING FUNCTION

build_hmms() {
    local generated_hmms=()
    for aln in "${aln_files[@]}"; do
        base_name=$(basename "$aln")
        prefix="${base_name%.*}"
        hmm_out="$outdir/${prefix}.hmm"
        log_file="$outdir/${prefix}_${timestamp}.log"

        echo "[INFO] Processing $base_name..."

        # Analyze alignment and get recommendations
        stats_tmp=$(mktemp)
        auto_params=($(analyze_alignment "$aln" 2> "$stats_tmp"))
        stats_info=$(cat "$stats_tmp")
        rm "$stats_tmp"

        cmd=(hmmbuild --cpu "$cpu")

        [[ "$wblosum" == true ]] && cmd+=(--wblosum)
        [[ "$wpb" == true ]] && cmd+=(--wpb)
        [[ -n "$wid" ]] && cmd+=(--wid "$wid")

        cmd+=("$hmm_out" "$aln")

        # 1. Logic to determine symfrac
        local sym_desc="Default (0.5)"
        if [[ -n "$symfrac" ]]; then
            cmd+=(--symfrac "$symfrac")
            sym_desc="Manual ($symfrac)"
        else
            for i in "${!auto_params[@]}"; do
                if [[ "${auto_params[$i]}" == "--symfrac" ]]; then
                    cmd+=(--symfrac "${auto_params[$((i+1))]}")
                    sym_desc="Recommended (${auto_params[$((i+1))]})"
                fi
            done
        fi

        # 2. Logic to determine weighting method
        local weight_desc="Default (GSC)"
        if [[ "$wblosum" == true ]]; then 
            cmd+=(--wblosum); weight_desc="Manual (--wblosum)" # If user explicitly chose wblosum 
        elif [[ "$wpb" == true ]]; then 
            cmd+=(--wpb); weight_desc="Manual (--wpb)" # If user explicitly chose wpb
        else
            for p in "${auto_params[@]}"; do
                if [[ "$p" == "--wblosum" ]]; then
                    cmd+=(--wblosum); weight_desc="Recommended (--wblosum)" # If infoalign recommended wblosum based on low identity
                elif [[ "$p" == "--wpb" ]]; then
                    cmd+=(--wpb); weight_desc="Recommended (--wpb)" # If infoalign recommended wpb based on high identity
                fi
            done
        fi

        [[ "$hand" == true ]] && cmd+=(--hand)
        cmd+=("$hmm_out" "$aln") # hmmbuild syntax: hmmbuild [options] <hmmfile> <alnfile>

        if [[ "$dry_run" == true ]]; then
            echo -e "\n\033[1;33m=== Recomendation(DRY-RUN) ===\033[0m"
            echo -e "Analysis Stats: $stats_info"
            echo -e "Schematic Weighting: $weight_desc"
            echo -e "Symfrac: $sym_desc"
            echo -e "Command that would be executed: \033[1;32m${cmd[*]}\033[0m\n"
        else
            {
                echo "=== HMMER PIPELINE LOG ==="
                echo "Date: $(date)"
                echo "Analysis Stats: $stats_info"
                echo "Weighting Scheme: $weight_desc"
                echo "Symfrac: $sym_desc"
                echo "Final Command: ${cmd[*]}"
                echo "---------------------------"
            } > "$log_file"

            echo -e "\033[0;32m[EXEC]\033[0m Generando perfil..."
            "${cmd[@]}" >> "$log_file" 2>&1
            generated_hmms+=("$hmm_out")
        fi
    done
    hmm_files+=("${generated_hmms[@]}")
}

# HMMSEARCH FUNCTION

run_search() {
    echo "[INFO] Running HMM search..."
    for hmm in "${hmm_files[@]}"; do
        prefix=$(basename "$hmm" .hmm)
        # search file log created in the back step, but we can also create a search-specific log if needed
        log_file="$outdir/${prefix}_${timestamp}.log"

        # hmm_base=$(basename "$hmm")
        # hmm_id="${hmm_base//./_}"
        # hmm_id="${hmm_id%.hmm}"

        for db in "${db_files[@]}"; do
            db_prefix=$(basename "$db" | cut -d. -f1)
            {
                echo -e "\n=== SEARCH STEP ==="
                echo "Database: $db"
                echo "Time: $(date)"
            } >> "$log_file"

            # db_base=$(basename "$db")
            # db_id="${db_base//./_}"
            # db_id="${db_id%.*}"

            cmd=()
            [[ "$db_type" == "prot" ]] && cmd=(hmmsearch --cpu "$cpu") || cmd=(nhmmer --cpu "$cpu")
            

            # if [[ "$db_type" == "prot" ]]; then
            #     cmd=(hmmsearch --cpu "$cpu")
            # else
            #     cmd=(nhmmer --cpu "$cpu")
            # fi

            [[ -n "$evalue" ]] && cmd+=(--E "$evalue")
            [[ -n "$domevalue" ]] && cmd+=(--domE "$domevalue")
            [[ "$tblout" == true ]] && cmd+=(--tblout "$outdir/${prefix}_vs_${db_prefix}.tbl")
            [[ "$domtblout" == true ]] && cmd+=(--domtblout "$outdir/${prefix}_vs_${db_prefix}.domtbl")
            [[ "$pfamtblout" == true ]] && cmd+=(--pfamtblout "$outdir/${prefix}_vs_${db_prefix}_pfamtblout.txt")

            cmd+=("$hmm" "$db")
            
            echo "Search Command: ${cmd[*]}" >> "$log_file"
            "${cmd[@]}" >> "$log_file" 2>&1
        done
    done
}

# EXECUTION STEPS

[[ "$mode" == "build" || "$mode" == "all" ]] && build_hmms
[[ "$mode" == "search" || "$mode" == "all" ]] && run_search

echo "[FINISHED] Logs y resultados en: $outdir"
# echo "[ALL DONE] Results saved in '$outdir'"
# echo "[SUMMARY]"
# echo "  HMMs built: ${#hmm_files[@]}"
# echo "  Databases searched: ${#db_files[@]}"
# echo "  Output directories:"
# echo "    HMMs:    $hmms_dir"
# echo "    Search:  $search_dir"
# echo "    Logs:    $logs_dir"

