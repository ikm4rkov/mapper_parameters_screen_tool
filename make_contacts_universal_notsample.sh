#!/usr/bin/env bash

# CONTACT_SCRIPT="/gpfs/ryabykhgrigory/bin/Bam_to_Contacts_grisha_fix.py"
# BASE_DIR="bowtie2_montecarlo_wrapper_outputs"
# MODE="BWA"
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
    echo "  -k, --mode HISAT|STAR|BWA (BWA also for bowtie2)"
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

# =========================================================
# Argument parsing
# =========================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--contact-mode)  MODE_SELECTION="$2"; shift 2 ;;
        -k)                 RNA_MODE="$2"; shift 2 ;;
        -kk)                DNA_MODE="$2"; shift 2 ;;
        --contact_script)   CONTACT_SCRIPT="$2"; shift 2 ;;
        --base-rna-dir)      BASE_RNA_DIR="$2"; shift 2 ;;
        --base-dna-dir)      BASE_DNA_DIR="$2"; shift 2 ;;
        --raw-contacts-dir) RAW_CONTACTS_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$MODE_SELECTION" ]]; then
    echo "Error: --mode is required."
    usage
fi

# =========================================================
# Collect all files
# =========================================================
declare -a DNA_FILES
declare -a RNA_FILES
declare -a MATCHED_PAIRS

# Собираем все файлы ДНК и РНК
mapfile -t DNA_FILES < <(find "$BASE_DNA_DIR" -type f -path "*dna_*/*_sorted.bam" | sort)
mapfile -t RNA_FILES < <(find "$BASE_RNA_DIR" -type f -path "*rna_*/*_sorted.bam" | sort)

# Строим пары по совпадению имени файла
# Строим все возможные пары ДНК с каждой РНК
for dna_file in "${DNA_FILES[@]}"; do
    for rna_file in "${RNA_FILES[@]}"; do
        MATCHED_PAIRS+=("${dna_file}:${rna_file}")
    done
done

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
fi

if [[ "${RNA_MODE^^}" == "BOWTIE2" ]]; then
    RNA_MODE="STAR" #HISAT
fi
if [[ "${DNA_MODE^^}" == "BOWTIE2" ]]; then
    DNA_MODE="STAR" #HISAT
fi
if [[ "${RNA_MODE^^}" == "HISAT2" ]]; then
    RNA_MODE="HISAT"
fi
if [[ "${DNA_MODE^^}" == "HISAT2" ]]; then
    DNA_MODE="HISAT"
fi
# =========================================================
# Run contact generation
# =========================================================
echo "Running contact generation on ${#SELECTED[@]} pairs..."
echo "Mapper-branch: $MAPPER_BRANCH"

    # Проверяем и создаём RAW_CONTACTS_DIR, если нужно
    if [[ -n "$RAW_CONTACTS_DIR" ]]; then
        mkdir -p "$RAW_CONTACTS_DIR"
    else
        echo "Error: --raw-contacts-dir must be specified."
        exit 1
    fi
    
    for entry in "${SELECTED[@]}"; do
        IFS=":" read -r dna_file rna_file <<< "$entry"
        
        # Формируем уникальный префикс для каждой пары
        dna_dir=$(basename "$(dirname "$dna_file")")
        rna_dir=$(basename "$(dirname "$rna_file")")
        p_arg="${rna_dir}__${dna_dir}"
        
        out_dir="$RAW_CONTACTS_DIR/$p_arg"
        cd "$RAW_CONTACTS_DIR"
        mkdir -p "$out_dir"
        cd "$out_dir"
        
        echo "DNA: $dna_file"
        echo "RNA: $rna_file"
        echo "Output prefix: $p_arg"
        
        python "$CONTACT_SCRIPT" \
            -r1 "$rna_file" \
            -r2 "$dna_file" \
            -mr "$RNA_MODE" \
            -md "$DNA_MODE" \
            -e "$ENZYME" \
            -p "raw_contacts"
    done
    
echo "All processing complete."

