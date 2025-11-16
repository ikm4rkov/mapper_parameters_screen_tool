#!/bin/bash
# Run EditDistance_CIGAR_filter.py on all *.tab.rc files in the given directory

FILTER_SCRIPT="/gpfs/ryabykhgrigory/bin/EditDistance_CIGAR_filter.py"
INPUT_DIR="/gpfs/ryabykhgrigory/mapping_hyperparameters_scan/bowtie_40outof50_rawcontacts"
OUTPUT_BASE="$INPUT_DIR/filtered"

# Parameters
MODE="NM + N_softClipp_bp"
PARAMS=(1 1 0 0 200 "no" "explorer" "ATA, not iMARGI")

mkdir -p "$OUTPUT_BASE"

if [[ ! -f "$FILTER_SCRIPT" ]]; then
    echo "❌ Error: filter script not found at $FILTER_SCRIPT"
    exit 1
fi

cd "$INPUT_DIR" || { echo "❌ Cannot access $INPUT_DIR"; exit 1; }

for file in "$INPUT_DIR"/*.tab.rc; do
    [[ -e "$file" ]] || { echo "No .tab.rc files found in $INPUT_DIR"; exit 1; }

    file_name=$(basename "$file")
    dir_name="${file_name%.tab.rc}"
    out_dir="$OUTPUT_BASE/$dir_name"
    mkdir -p "$out_dir"

    echo "🔹 Processing $file_name ..."
    echo "Command: python3 $FILTER_SCRIPT \"$MODE\" ${PARAMS[*]} \"$file_name\" \"$INPUT_DIR/\" \"$out_dir/\""
    
    python3 "$FILTER_SCRIPT" "$MODE" "${PARAMS[@]}" "$file_name" "$INPUT_DIR/" "$out_dir/"

    if compgen -G "$out_dir/*filtered*.tsv" > /dev/null; then
        echo "✅ Done: output saved in $out_dir/"
    else
        echo "⚠️  Warning: no filtered file detected for $file_name"
    fi
done

echo "All files processed."

