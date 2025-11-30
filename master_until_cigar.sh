#!/usr/bin/env bash
set -euo pipefail

###############################################
# Default configuration
###############################################

BASE_DIR=""
SCRIPT2_MODE_SELECTION="all"
SCRIPT2_COUNT=0
SCRIPT2_SEED=0
SCRIPT2_DIFFERENT=0
SCRIPT2_MAPPER_BRANCH="collab"

SAMPLE_CONTACTS=0
CONTACTS_MODE="all"
CONTACTS_PAIRS_FILE="pairs.tsv"

RAW_CONTACTS_DIR=""
FILTER_INPUT_DIR=""
FILTER_OUTPUT_DIR=""
FILTER_RESULTS_FILE="results_summary.txt"

LOG_FILE="wrapper.log"

###############################################
# Usage function
###############################################
usage() {
    echo "Usage: $0 -b <base_dir> [options]"
    echo
    echo "Options:"
    echo "  -b, --base-dir DIR          Base directory for BAM/SAM processing (required)"
    echo "  --sample-contacts 0|1       Whether to run sample contacts script (default 0)"
    echo "  --mode all|one              Mode for contacts generation script (default all)"
    echo "  --pairs FILE                Pairs TSV file for contacts script (default pairs.tsv)"
    echo "  --raw-contacts DIR          Directory to save raw contact files"
    echo "  --mapper-branch native|collab  Branch for script #2 (default collab)"
    echo
    echo "Options for script #2 random mode:"
    echo "  -c, --count NUMBER"
    echo "  -s, --seed NUMBER"
    echo "  -d, --different 0|1"
    exit 1
}

###############################################
# Parse command-line arguments
###############################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--base-dir) BASE_DIR="$2"; shift 2 ;;
        --sample-contacts) SAMPLE_CONTACTS="$2"; shift 2 ;;
        --mode) CONTACTS_MODE="$2"; shift 2 ;;
        --pairs) CONTACTS_PAIRS_FILE="$2"; shift 2 ;;
        --raw-contacts) RAW_CONTACTS_DIR="$2"; shift 2 ;;
        --mapper-branch) SCRIPT2_MAPPER_BRANCH="$2"; shift 2 ;;
        -c|--count) SCRIPT2_COUNT="$2"; shift 2 ;;
        -s|--seed) SCRIPT2_SEED="$2"; shift 2 ;;
        -d|--different) SCRIPT2_DIFFERENT="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$BASE_DIR" ]]; then
    echo "Error: --base-dir is required."
    usage
fi

###############################################
# Validate raw contacts directory
###############################################
if [[ -n "$RAW_CONTACTS_DIR" ]]; then
    if [[ -e "$RAW_CONTACTS_DIR" && ! -d "$RAW_CONTACTS_DIR" ]]; then
        echo "Error: --raw-contacts-dir exists but is not a directory: $RAW_CONTACTS_DIR" >&2
        exit 1
    fi
    mkdir -p "$RAW_CONTACTS_DIR"
    if [[ ! -w "$RAW_CONTACTS_DIR" ]]; then
        echo "Error: Directory not writable: $RAW_CONTACTS_DIR" >&2
        exit 1
    fi
fi

echo "Logging to $LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

###############################################
# Step 1: Sort BAM/SAM files
###############################################
./sam2bam_and_sort.sh "$BASE_DIR"

###############################################
# Step 2: Generate contacts
###############################################
if [[ "$SAMPLE_CONTACTS" -eq 1 ]]; then
    ./make_contacts_sample.sh --mode "$CONTACTS_MODE" --pairs "$CONTACTS_PAIRS_FILE"
else
    SCRIPT2_ARGS=(--base-dir "$BASE_DIR" -m "$CONTACTS_MODE" --mapper-branch "$SCRIPT2_MAPPER_BRANCH")
    if [[ "$SCRIPT2_COUNT" -gt 0 ]]; then
        SCRIPT2_ARGS+=(-m random -c "$SCRIPT2_COUNT" -s "$SCRIPT2_SEED" -d "$SCRIPT2_DIFFERENT")
    fi
    ./make_contacts_universal_notsample.sh "${SCRIPT2_ARGS[@]}"
fi

###############################################
# Step 2b: Move raw contacts
###############################################
if [[ -n "$RAW_CONTACTS_DIR" ]]; then
    mv ./*.tab.rc "$RAW_CONTACTS_DIR"/ 2>/dev/null || true
fi

###############################################
# Step 3: Filter contacts
###############################################
if [[ -n "$RAW_CONTACTS_DIR" ]]; then
    FILTER_INPUT_DIR="$RAW_CONTACTS_DIR"
fi

FILTER_OUTPUT_DIR="${FILTER_OUTPUT_DIR:-filtered_contacts}"
mkdir -p "$FILTER_OUTPUT_DIR"

if [[ -n "$FILTER_INPUT_DIR" ]]; then
    ./calculate_CIAGAR_filter.sh -i "$FILTER_INPUT_DIR" -o "$FILTER_OUTPUT_DIR" -r "$FILTER_RESULTS_FILE"
fi


