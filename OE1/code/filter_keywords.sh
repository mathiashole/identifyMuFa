#!/bin/bash

# Function to check if a file exists and is readable
check_file_exists() {
    local file="$1"

    if [ ! -f "$file" ]; then
        echo "Error: File \"$file\" not found!"
        return 1
    fi

    if [ ! -r "$file" ]; then
        echo "Error: File \"$file\" is not readable!"
        return 1
    fi

    return 0
}

# Main script starts here

# Check if input arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_table> <output_dir>"
    exit 1
fi

input_file="$1"
output_dir="$2"

# Check if input file exists and is readable
if ! check_file_exists "$input_file"; then
    exit 1
fi

# Create the output directory if it doesn't exist
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    echo "Created output directory: $output_dir"
fi

# Print content of input file
echo "Content of $input_file:"
cat "$input_file"
echo "-------------------------"

# Read each line from input file and process
while IFS=$'\t' read -r _ file keyword1 keyword2 || [ -n "$file" ]; do
    # Remove extra double quotes from file path
    file="${file//\"/}"

    # Remove double quotes from keywords
    keyword1="${keyword1//\"/}"
    keyword2="${keyword2//\"/}"

    # Check if the file exists and is readable
    if ! check_file_exists "$file"; then
        continue
    fi

    # Extract filename without directory
    filename=$(basename "$file")

    # Generate output file name in the output directory
    output_file="$output_dir/filtered_:${keyword1}_${keyword2}:_$filename"

    # Generate output no filter file name in the output directory
    output_file="$output_dir/non-filtered_:${keyword2}:_$filename"

    # Debugging output
    echo "Processing file: $file"
    echo "Keywords: $keyword1, $keyword2"

    # Use grep to filter the content based on keywords
    grep "$keyword1" "$file" | grep "$keyword2" | sort | uniq > "$output_file"

    # Check if output file is not empty
    if [ -s "$output_file" ]; then
        echo "Filtered data saved to \"$output_file\""
    else
        echo "No data found for \"$keyword1\" and \"$keyword2\" in \"$file\""
    fi
done < "$input_file"
