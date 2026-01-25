#!/bin/bash

# ---------------------------
# Parse command-line options
# ---------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) input_dir="$2"; shift 2 ;;
        -o) output_base="$2"; shift 2 ;;
        -s) EDIT_DISTANCE_CIGAR_FILTER_SCRIPT="$2"; shift 2 ;;
        -k)                 RNA_MODE="$2"; shift 2 ;;
        -kk)                DNA_MODE="$2"; shift 2 ;;
        *)
            echo "Usage: $0 -i input_dir -o output_base -s EDIT_DISTANCE_CIGAR_FILTER_SCRIPT"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$input_dir" ] || [ -z "$output_base" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 -i input_dir -o output_base"
    exit 1
fi

# ---------------------------
# Main script logic
# ---------------------------

# Ensure output directories exist
mkdir -p "$output_base"

# Process all *.tab.rc files
for subdir in "$input_dir"/rna_${RNA_MODE,,}_*__dna_${DNA_MODE,,}_*/; do
    # Получаем имя поддиректории
    subdir_name=$(basename "$subdir")
    # Создаём соответствующую папку в output_base
    output_dir="$output_base/$subdir_name"
    mkdir -p "$output_dir"
    echo "Processing directory: $subdir_name"
    # # Full path to input file
    # input_file_path="$input_dir/$file_name"

    # Run Python script
    python "$EDIT_DISTANCE_CIGAR_FILTER_SCRIPT" \
        "NM + N_softClipp_bp" 2 2 0 0 200 "no" "explorer" "ATA, not iMARGI" \
        "raw_contacts_Unique_RNA.tab.rc" \
        "$subdir/" \
        "$output_dir/"

done

echo "CIGAR filter complete."