#!/bin/bash

# ---------------------------
# Parse command-line options
# ---------------------------
while getopts "i:o:r:" opt; do
    case "$opt" in
        i) input_dir="$OPTARG" ;;
        o) output_base="$OPTARG" ;;
        r) results_file="$OPTARG" ;;
        *)
            echo "Usage: $0 -i input_dir -o output_base -r results_file"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$input_dir" ] || [ -z "$output_base" ] || [ -z "$results_file" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 -i input_dir -o output_base -r results_file"
    exit 1
fi

# ---------------------------
# Main script logic
# ---------------------------

# Ensure output directories exist
mkdir -p "$output_base"

# Clear or create results file
> "$results_file"

# Move into output base directory
cd "$output_base" || exit 1

# Process all *.tab.rc files
for input_file in "$input_dir"/*Unique_RNA.tab.rc; do
    file_name=$(basename "$input_file")

    # Create a subdirectory for this file
    output_dir="$output_base/${file_name%.tab.rc}"
    mkdir -p "$output_dir"

    # Full path to input file
    input_file_path="$input_dir/$file_name"

    # Run Python script
    python /gpfs/ryabykhgrigory/bin/EditDistance_CIGAR_filter.py \
        "NM + N_softClipp_bp" 2 2 0 0 200 "no" "explorer" "ATA, iMARGI" \
        "$file_name" \
        "${input_dir}/" \
        "$output_dir"

    # Log filename to results
    echo "$file_name" >> "$results_file"

    # Path to filtered file
    filtered_file="$output_dir/filtered_${file_name}"

    # Analyze if exists
    if [ -f "$filtered_file" ]; then
        cut -f2 "$filtered_file" | sort | uniq -c >> "$results_file"
    else
        echo "File $filtered_file not found" >> "$results_file"
    fi

    echo "" >> "$results_file"
done

echo "Analysis complete. Results saved to $results_file"

