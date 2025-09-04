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
        avg_id = (n > 0 ? sum_id/n : 0)
        avg_len = (n > 0 ? sum_len/n : 0)

        for(i=1;i<=n;i++){
            sd_id += (id[i]-avg_id)^2
            sd_len += (len[i]-avg_len)^2
        }
        sd_id = (n>1 ? sqrt(sd_id/(n-1)) : 0)
        sd_len = (n>1 ? sqrt(sd_len/(n-1)) : 0)

        asort(id)
        asort(len)
        if(n % 2){
            med_id=id[(n+1)/2]
            med_len=len[(n+1)/2]
        } else {
            med_id=(id[n/2]+id[n/2+1])/2
            med_len=(len[n/2]+len[n/2+1])/2
        }

        print "Queries total:", length(queries)
        print "Queries with hits:", length(qh)
        print "Unique queries:", length(queries)
        print "Unique subjects:", length(subjects)
        print "Avg identity:", avg_id
        print "SD identity:", sd_id
        print "Median identity:", med_id
        print "Avg length:", avg_len
        print "SD length:", sd_len
        print "Median length:", med_len
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
    }' "$BLAST_FILE" | sort -k1,1 -k3,3nr -k4,4nr | "$BLAST_FILE" | awk '!seen[$1]++' > "$OUTPUT_FILE"

else
    awk -v id_min="$IDENTITY_MIN" -v len_min="$LENGTH_MIN" -v evalue_max="$EVALUE_MAX" -v col="$COLUMN_TO_PRINT" '
    BEGIN { FS=OFS="\t" }
    ($3 >= id_min) && ($4 >= len_min) && ($11 <= evalue_max) && (self_filter == "false" || $1 != $2) {
        if (col == "ALL") print $0;
        else print $col;
    }' "$BLAST_FILE" > "$OUTPUT_FILE"
fi


if [[ $REPORT == true ]]; then
    generate_report "$INPUT_FILE" "Raw data (unfiltered)"
    generate_report "$OUTPUT_FILE" "Filtered data"
fi

echo "Blast Refiner results saved to: $OUTPUT_FILE"
