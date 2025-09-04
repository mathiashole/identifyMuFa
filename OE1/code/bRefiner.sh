#!/bin/bash

# Default values
IDENTITY_MIN=0
LENGTH_MIN=0
EVALUE_MAX=1  # Default 1 (accepts all values)
UNIQ_SORT=false    # By default, do NOT apply uniq
BEST_ALIGNMENT=false  # By default, do NOT pick best hit
COLUMN_TO_PRINT="ALL" # Default: print all columns
SELF_FILTER=true   # By default, filter self-hits
REPORT=false # By default report
Q_COV_MIN=0
S_COV_MIN=0

generate_report() {
    local file="$1"
    awk '
    BEGIN { FS="\t"; OFS="\t"; }
    NR==1 {next} # skip header if any
    {
        q[$1]++; s[$2]++;
        id_sum+=$3; id2_sum+=$3*$3;
        len_sum+=$4; len2_sum+=$4*$4;
        n++
    }
    END {
        if (n==0) {print "No data"; exit}
        avg_id=id_sum/n
        avg_len=len_sum/n
        var_id=(id2_sum/n - avg_id*avg_id)
        var_len=(len2_sum/n - avg_len*avg_len)
        sd_id=sqrt(var_id)
        sd_len=sqrt(var_len)
        print "Queries total:", n
        print "Unique queries:", length(q)
        print "Unique subjects:", length(s)
        print "Avg identity:", avg_id
        print "Std identity:", sd_id
        print "Avg length:", avg_len
        print "Std length:", sd_len
    }' "$file"
    echo "------------------------"
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -file) BLAST_FILE="$2"; shift ;;
        -i) IDENTITY_MIN="$2"; shift ;;
        -l) LENGTH_MIN="$2"; shift ;;
        -e) EVALUE_MAX="$2"; shift ;;  # Add evalue filter
        -col) COLUMN_TO_PRINT="$2"; shift ;;
        -uniq) UNIQ_SORT=true ;;  # Apply sorting and unique filtering
        -best) BEST_ALIGNMENT=true ;;
        -self) SELF_FILTER="$2"; shift ;;  # Add self-hits filter option
        -report) REPORT=true ;; # Apply report option
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$BLAST_FILE" ]] || [[ -z "$LENGTH_MIN" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -file <filename> -l <length_min> [-e <evalue_max>] [-col <column_to_print>] [-i <identity_min>] [-uniq] [-best] [-self <true|false>]"
    exit 1
fi

# Get the path and base name of the file
INPUT_DIR=$(dirname "$BLAST_FILE")
INPUT_NAME=$(basename "$BLAST_FILE")

# Construct the output file name
OUTPUT_FILE="$INPUT_DIR/bRefiner_$INPUT_NAME"

# Apply filters
if [ "$UNIQ_SORT" = true ] && [ "$BEST_ALIGNMENT" = true ]; then
    echo "Error: -uniq and -best cannot be used together."
    exit 1
fi

# Apply filters
if [ "$UNIQ_SORT" = true ]; then
    awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v col="$COLUMN_TO_PRINT" '
    BEGIN { FS=OFS="\t" }
    ($3 >= id_min) && ($4 >= len_min) && ($11 <= evalue_max) && (self_filter == "false" || $1 != $2) {
        if (col == "ALL") print $0;
        else print $col;
    }' "$BLAST_FILE" | sort | uniq > "$OUTPUT_FILE"

elif [ "$BEST_ALIGNMENT" != "" ]; then
    awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v col="$COLUMN_TO_PRINT" '
    BEGIN { FS=OFS="\t" }
    ($3 >= id_min) && ($4 >= len_min) && ($11 <= evalue_max) && (self_filter == "false" || $1 != $2) {
        if (col == "ALL") print $0;
        else print $col;
    }' "$BLAST_FILE" \
    | sort -k1,1 -k3,3nr -k4,4nr \
    | awk '!seen[$1]++' > "$OUTPUT_FILE"

else
    awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v col="$COLUMN_TO_PRINT" '
    BEGIN { FS=OFS="\t" }
    ($3 >= id_min) && ($4 >= len_min) && ($11 <= evalue_max) && (self_filter == "false" || $1 != $2) {
        if (col == "ALL") print $0;
        else print $col;
    }' "$BLAST_FILE" > "$OUTPUT_FILE"
fi

# Run report if requested
if [ "$REPORT" = true ]; then
    echo "Report for: Raw data (unfiltered)"
    generate_report "$BLAST_FILE"

    echo
    echo "Report for: Filtered data"
    generate_report "$OUTPUT_FILE"
fi
# if [[ $REPORT == true ]]; then
#     generate_report "$INPUT_FILE" "Raw data (unfiltered)"
#     generate_report "$OUTPUT_FILE" "Filtered data"
# fi

echo "Blast Refiner results saved to: $OUTPUT_FILE"
