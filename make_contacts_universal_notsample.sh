#!/usr/bin/env bash

CONTACT_SCRIPT="/gpfs/ryabykhgrigory/bin/Bam_to_Contacts_grisha_fix.py"
BASE_DIR="bowtie2_montecarlo_wrapper_outputs"
MODE="STAR"
ENZYME="ATA"

# =========================================================
# Usage
# =========================================================
usage() {
    echo "Usage: $0 -m <mode> [options]"
    echo
    echo "Modes:"
    echo "  all      Process all DNA/RNA pairs (default)"
    echo "  random   Process random subset"
    echo
    echo "Options for random mode:"
    echo "  -c, --count <number>     How many random pairs"
    echo "  -s, --seed <number>      Random seed"
    echo "  -d, --different <0|1>    0: same matched pairs, 1: random DNA/RNA separately"
    echo
    echo "General options:"
    echo "  -b, --base-dir <dir>     Override input directory"
    echo "      --mapper-branch <collab|native>"
    echo "                           collab = use folder names in prefix"
    echo "                           native = use BAM base filenames in prefix"
    echo
    exit 1
}

# =========================================================
# Defaults
# =========================================================
MODE_SELECTION=""
COUNT=0
SEED=0
DIFFERENT=0
BASE_DIR_OVERRIDE=""
MAPPER_BRANCH="collab"   # default behavior preserved

# =========================================================
# Argument parsing
# =========================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mode) MODE_SELECTION="$2"; shift 2 ;;
        -c|--count) COUNT="$2"; shift 2 ;;
        -s|--seed) SEED="$2"; shift 2 ;;
        -d|--different) DIFFERENT="$2"; shift 2 ;;
        -b|--base-dir) BASE_DIR_OVERRIDE="$2"; shift 2 ;;
        --mapper-branch)
            MAPPER_BRANCH="$2"
            if [[ "$MAPPER_BRANCH" != "collab" && "$MAPPER_BRANCH" != "native" ]]; then
                echo "Error: --mapper-branch must be 'collab' or 'native'"
                exit 1
            fi
            shift 2
            ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$MODE_SELECTION" ]]; then
    echo "Error: --mode is required."
    usage
fi

# =========================================================
# Override BASE_DIR if requested
# =========================================================
if [[ -n "$BASE_DIR_OVERRIDE" ]]; then
    if [[ ! -d "$BASE_DIR_OVERRIDE" ]]; then
        echo "Error: Provided --base-dir '$BASE_DIR_OVERRIDE' is not a directory."
        exit 1
    fi
    BASE_DIR="$BASE_DIR_OVERRIDE"
fi

if [[ "$MODE_SELECTION" == "random" && ( "$COUNT" -le 0 || "$SEED" -eq 0 ) ]]; then
    echo "Error: --count and --seed must be positive for random mode."
    usage
fi

if [[ "$MODE_SELECTION" == "random" && "$DIFFERENT" -ne 0 && "$DIFFERENT" -ne 1 ]]; then
    echo "Error: --different must be 0 or 1."
    usage
fi

# =========================================================
# Deterministic random generator
# =========================================================
getseededrandom() {
    local seed="$1"
    awk -v seed="$seed" 'BEGIN { srand(seed); for(i=0;i<65536;i++) printf("%c", int(rand()*256)); }'
}

# =========================================================
# Collect all files
# =========================================================
declare -a DNA_FILES
declare -a RNA_FILES
declare -a MATCHED_PAIRS

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
# Select pairs
# =========================================================
if [[ "$MODE_SELECTION" == "all" ]]; then
    SELECTED=("${MATCHED_PAIRS[@]}")

elif [[ "$MODE_SELECTION" == "random" && "$DIFFERENT" -eq 0 ]]; then
    if (( COUNT > total_pairs )); then
        echo "Requested $COUNT pairs but only $total_pairs available."
        exit 1
    fi

    mapfile -t SELECTED < <(
        printf "%s\n" "${MATCHED_PAIRS[@]}" |
        shuf --random-source=<(getseededrandom "$SEED") |
        head -n "$COUNT"
    )

elif [[ "$MODE_SELECTION" == "random" && "$DIFFERENT" -eq 1 ]]; then
    if (( COUNT > total_dna || COUNT > total_rna )); then
        echo "Requested $COUNT pairs but only $total_dna DNA and $total_rna RNA files available."
        exit 1
    fi

    mapfile -t DNA_SEL < <(
        printf "%s\n" "${DNA_FILES[@]}" |
        shuf --random-source=<(getseededrandom "$SEED") |
        head -n "$COUNT"
    )

    mapfile -t RNA_SEL < <(
        printf "%s\n" "${RNA_FILES[@]}" |
        shuf --random-source=<(getseededrandom "$((SEED+1))") |
        head -n "$COUNT"
    )

    SELECTED=()
    for idx in $(seq 0 $((COUNT-1))); do
        SELECTED+=("${DNA_SEL[$idx]}:${RNA_SEL[$idx]}")
    done
fi

# =========================================================
# Run contact generation
# =========================================================
echo "Running contact generation on ${#SELECTED[@]} pairs..."
echo "Mapper-branch: $MAPPER_BRANCH"

for entry in "${SELECTED[@]}"; do
    IFS=":" read -r dna_file rna_file <<< "$entry"

    # --- collab mode (original behavior) ---
    if [[ "$MAPPER_BRANCH" == "collab" ]]; then
        dna_dir=$(basename "$(dirname "$dna_file")")
        rna_dir=$(basename "$(dirname "$rna_file")")
        base_name=$(basename "$dna_file" "_sorted.bam")
        p_arg="${dna_dir}__${rna_dir}__${base_name}"

    # --- native mode (new behavior) ---
    else
        dna_base=$(basename "$dna_file" ".bam")
        rna_base=$(basename "$rna_file" ".bam")
        p_arg="${dna_base}__${rna_base}"
    fi

    echo
    echo "DNA: $dna_file"
    echo "RNA: $rna_file"
    echo "Output prefix: $p_arg"

    python "$CONTACT_SCRIPT" \
        -r1 "$rna_file" \
        -r2 "$dna_file" \
        -mr "$MODE" \
        -md "$MODE" \
        -e "$ENZYME" \
        -p "$p_arg"
done

echo "All processing complete."

