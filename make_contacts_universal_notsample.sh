#!/usr/bin/env bash

CONTACT_SCRIPT="/gpfs/ryabykhgrigory/bin/Bam_to_Contacts_grisha_fix.py"
BASE_DIR="hisat2_test_output"
MODE="STAR"
ENZYME="ATA"

# =========================================================
# Argument parsing
# =========================================================
usage() {
    echo "Usage: $0 -m <mode> [options]"
    echo
    echo "Modes:"
    echo "  all      Process all DNA/RNA pairs (default behavior)."
    echo "  random   Process random subset of pairs."
    echo
    echo "Options for random mode:"
    echo "  -c, --count <number>     Number of random pairs to process."
    echo "  -s, --seed <number>      Random seed for reproducibility."
    echo "  -d, --different <0|1>    0 = same pairs, 1 = random DNA and RNA separately."
    echo
    echo "Examples:"
    echo "  $0 -m all"
    echo "  $0 -m random -c 250 -s 1"
    echo "  $0 -m random -c 250 -s 1 -d 1"
    exit 1
}

# Defaults
MODE_SELECTION=""
COUNT=0
SEED=0
DIFFERENT=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mode) MODE_SELECTION="$2"; shift 2 ;;
        -c|--count) COUNT="$2"; shift 2 ;;
        -s|--seed) SEED="$2"; shift 2 ;;
        -d|--different) DIFFERENT="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$MODE_SELECTION" ]]; then
    echo "Error: --mode is required."
    usage
fi

if [[ "$MODE_SELECTION" == "random" && ( "$COUNT" -le 0 || "$SEED" -eq 0 ) ]]; then
    echo "Error: --count and --seed must be specified and positive for random mode."
    usage
fi

if [[ "$MODE_SELECTION" == "random" && "$DIFFERENT" -ne 0 && "$DIFFERENT" -ne 1 ]]; then
    echo "Error: --different must be 0 or 1."
    usage
fi

# =========================================================
# Helper: deterministic random generator for shuf
# =========================================================
getseededrandom() {
    local seed="$1"
    awk -v seed="$seed" 'BEGIN { srand(seed); for(i=0;i<65536;i++) printf("%c", int(rand()*256)); }'
}

# =========================================================
# Collect all valid DNA/RNA files
# =========================================================
declare -a DNA_FILES
declare -a RNA_FILES
declare -a MATCHED_PAIRS

# Find all DNA directories (contain 'dna')
while IFS= read -r dna_dir; do
    [[ -d "$dna_dir" ]] || continue
    rna_dir="${dna_dir/dna/rna}"

    [[ -d "$rna_dir" ]] || continue

    while IFS= read -r dna_file; do
        [[ "$dna_file" == *"_sorted.bam" ]] || continue
        filename=$(basename "$dna_file")
        rna_file="${rna_dir}/${filename}"

        if [[ -f "$rna_file" ]]; then
            MATCHED_PAIRS+=("${dna_file}:${rna_file}")
        fi
        DNA_FILES+=("$dna_file")
    done < <(find "$dna_dir" -type f -name "*_sorted.bam" | sort)

    while IFS= read -r rna_file; do
        [[ "$rna_file" == *"_sorted.bam" ]] || continue
        RNA_FILES+=("$rna_file")
    done < <(find "$rna_dir" -type f -name "*_sorted.bam" | sort)
done < <(find "$BASE_DIR" -type d -name "*dna*" | sort)

total_pairs=${#MATCHED_PAIRS[@]}
total_dna=${#DNA_FILES[@]}
total_rna=${#RNA_FILES[@]}

echo "Found $total_dna DNA BAMs and $total_rna RNA BAMs."
echo "Detected $total_pairs valid matched pairs."

if [[ "$total_pairs" -eq 0 ]]; then
    echo "No valid DNA/RNA pairs found!"
    exit 1
fi

# =========================================================
# Mode selection
# =========================================================
if [[ "$MODE_SELECTION" == "all" ]]; then
    echo "Processing all $total_pairs pairs..."
    SELECTED=("${MATCHED_PAIRS[@]}")

elif [[ "$MODE_SELECTION" == "random" && "$DIFFERENT" -eq 0 ]]; then
    if (( COUNT > total_pairs )); then
        echo "Requested $COUNT pairs but only $total_pairs available."
        exit 1
    fi
    echo "Selecting $COUNT random *matched* pairs with seed $SEED..."
    mapfile -t SELECTED < <(printf "%s\n" "${MATCHED_PAIRS[@]}" |
        shuf --random-source=<(getseededrandom "$SEED") | head -n "$COUNT")

elif [[ "$MODE_SELECTION" == "random" && "$DIFFERENT" -eq 1 ]]; then
    if (( COUNT > total_dna || COUNT > total_rna )); then
        echo "Requested $COUNT pairs but only $total_dna DNA and $total_rna RNA files available."
        exit 1
    fi
    echo "Selecting $COUNT random DNA (seed=$SEED) and RNA (seed=$((SEED+1))) files..."
    mapfile -t DNA_SEL < <(printf "%s\n" "${DNA_FILES[@]}" |
        shuf --random-source=<(getseededrandom "$SEED") | head -n "$COUNT")
    mapfile -t RNA_SEL < <(printf "%s\n" "${RNA_FILES[@]}" |
        shuf --random-source=<(getseededrandom "$((SEED+1))") | head -n "$COUNT")

    SELECTED=()
    for idx in $(seq 0 $((COUNT-1))); do
        dna_file="${DNA_SEL[$idx]}"
        rna_file="${RNA_SEL[$idx]}"
        SELECTED+=("${dna_file}:${rna_file}")
    done
fi

# =========================================================
# Run contact generation
# =========================================================
echo "Running contact generation on ${#SELECTED[@]} selected entries..."

for entry in "${SELECTED[@]}"; do
    IFS=":" read -r dna_file rna_file <<< "$entry"

    dna_dir=$(basename "$(dirname "$dna_file")")
    rna_dir=$(basename "$(dirname "$rna_file")")
    base_name=$(basename "$dna_file" "_sorted.bam")
    p_arg="${dna_dir}__${rna_dir}__${base_name}"

    echo
    echo "Running contact generation for:"
    echo "  DNA: $dna_file"
    echo "  RNA: $rna_file"
    echo "  Output Prefix: $p_arg"

    python "$CONTACT_SCRIPT" \
        -r1 "$rna_file" \
        -r2 "$dna_file" \
        -mr "$MODE" \
        -md "$MODE" \
        -e "$ENZYME" \
        -p "$p_arg"
done

echo "All processing complete."

