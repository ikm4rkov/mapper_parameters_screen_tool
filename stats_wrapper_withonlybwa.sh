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
# Build job list (ONLY bwa-related pairs)
###############################################
echo "Scanning RNA/DNA directories (bwa-filtered)..."

JOBLIST=$(mktemp)

for rna_path in "$RAW_CONTACTS_DIR"/rna_*; do
    [ -d "$rna_path" ] || continue
    rna_dir=$(basename "$rna_path")

    for dna_path in "$RAW_CONTACTS_DIR"/dna_*; do
        [ -d "$dna_path" ] || continue
        dna_dir=$(basename "$dna_path")

        # Only include pairs where at least one contains "bwa"
        if [[ "$rna_dir" == *bwa* || "$dna_dir" == *bwa* ]]; then
            echo "$RAW_CONTACTS_DIR $FILTERED_CONTACTS_DIR $rna_dir $dna_dir" >> "$JOBLIST"
        fi
    done
done

###############################################
# Prepare tmp files (one per thread)
###############################################
echo "Preparing temporary files..."

for i in $(seq 1 "$THREADS"); do
    rm -f "$OUT_DIR/tmp_${i}.tsv"
done

###############################################
# Run in parallel
###############################################
echo "Running jobs with $THREADS threads..."

parallel -j "$THREADS" --colsep ' ' \
'
TMP_FILE="'$OUT_DIR'/tmp_{%}.tsv"

if [ ! -f "$TMP_FILE" ]; then
    # первый запуск в этом slot → пишем с header
    "'$EXE'" {1} {2} {3} {4} "$TMP_FILE"
else
    # последующие → убираем header
    "'$EXE'" {1} {2} {3} {4} /dev/stdout | tail -n +2 >> "$TMP_FILE"
fi
' :::: "$JOBLIST"

###############################################
# Merge outputs
###############################################
echo "Merging outputs..."

MERGED="$OUT_DIR/all_stats_merged.tsv"

TMP_FILES=("$OUT_DIR"/tmp_*.tsv)

if [[ ${#TMP_FILES[@]} -eq 0 ]]; then
    echo "No tmp files to merge."
    exit 1
fi

# первый файл с header
cat "${TMP_FILES[0]}" > "$MERGED"

# остальные без header
for f in "${TMP_FILES[@]:1}"; do
    tail -n +2 "$f" >> "$MERGED"
done

###############################################
# Cleanup
###############################################
rm -f "$OUT_DIR"/tmp_*.tsv
rm "$JOBLIST"

echo "All stats computed and merged into $MERGED."
