#!/usr/bin/env bash
set -euo pipefail

###############################################
# Usage
###############################################
usage() {
    echo "Usage:"
    echo "  $0 <config_file> [--exe contacts_stats] [-t threads]"
    exit 1
}

###############################################
# Parse config
###############################################
if [ $# -lt 1 ]; then
    usage
fi

config_file="$1"
shift

if [ ! -f "$config_file" ]; then
    echo "Error: Config file '$config_file' not found."
    exit 1
fi

# shellcheck source=/dev/null
source "$config_file"

###############################################
# Parse arguments
###############################################
EXE="./contacts_stats"
THREADS=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --exe) EXE="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

###############################################
# Directories
###############################################
RAW_CONTACTS_DIR="$WORK_DIR/raw_contacts"
FILTERED_CONTACTS_DIR="$WORK_DIR/filtered_contacts"
OUT_DIR="$WORK_DIR/intersection_stats"
mkdir -p "$OUT_DIR"

###############################################
# Build job list
###############################################
echo "Scanning RNA/DNA directories..."

JOBLIST=$(mktemp)

for rna_path in "$RAW_CONTACTS_DIR"/rna_*; do
    [ -d "$rna_path" ] || continue
    rna_dir=$(basename "$rna_path")

    for dna_path in "$RAW_CONTACTS_DIR"/dna_*; do
        [ -d "$dna_path" ] || continue
        dna_dir=$(basename "$dna_path")

        out_file="$OUT_DIR/stats_${rna_dir}__${dna_dir}.tsv"

        echo "$RAW_CONTACTS_DIR $FILTERED_CONTACTS_DIR $rna_dir $dna_dir $out_file" >> "$JOBLIST"
    done
done

###############################################
# Run in parallel
###############################################
echo "Running jobs with $THREADS threads..."

cat "$JOBLIST" | parallel -j "$THREADS" \
    "$EXE {1} {2} {3} {4} {5}"

###############################################
# Merge outputs
###############################################
echo "Merging outputs..."

MERGED="$OUT_DIR/all_stats_merged.tsv"

mapfile -t TMP_LIST < <(awk '{print $5}' "$JOBLIST")

if [[ ${#TMP_LIST[@]} -eq 0 ]]; then
    echo "No output files to merge."
    exit 1
fi

cat "${TMP_LIST[0]}" > "$MERGED"

for f in "${TMP_LIST[@]:1}"; do
    tail -n +2 "$f" >> "$MERGED"
done

rm "$JOBLIST"

echo "All stats computed and merged into $MERGED."
