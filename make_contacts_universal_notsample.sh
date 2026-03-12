#!/usr/bin/env bash

# =========================================================
# Usage
# =========================================================
# usage() {
#     echo "Usage: $0 -m <mode> [options]"
#     echo
#     echo "Modes:"
#     echo "  all      Process all DNA/RNA pairs (default)"
#     echo "  random   Process random subset"
#     echo
#     echo "Options for random mode:"
#     echo "  -c, --count <number>     How many random pairs"
#     echo "  -s, --seed <number>      Random seed"
#     echo "  -d, --different <0|1>    0: same matched pairs, 1: random DNA/RNA separately"
#     echo
#     echo "General options:"
#     echo "  -b, --base-dir <dir>     Override input directory"
#     echo "      --mapper-branch <collab|native>"
#     echo "                           collab = use folder names in prefix"
#     echo "                           native = use BAM base filenames in prefix"
#     echo "  -k, --mode HISAT|STAR|BWA (BWA also for bowtie2)"
#     echo
#     exit 1
# }


# =========================================================
# Argument parsing
# =========================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -k)                 RNA_MODE="$2"; shift 2 ;;
        -kk)                DNA_MODE="$2"; shift 2 ;;
        --contact_script)   CONTACT_SCRIPT="$2"; shift 2 ;;
        --base-rna-dir)      BASE_RNA_DIR="$2"; shift 2 ;;
        --base-dna-dir)      BASE_DNA_DIR="$2"; shift 2 ;;
        --raw-contacts-dir) RAW_CONTACTS_DIR="$2"; shift 2 ;;
        # -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# =========================================================
# Collect all files
# =========================================================
declare -a DNA_FILES
declare -a RNA_FILES

# Собираем все файлы ДНК и РНК
mapfile -t DNA_FILES < <(find "$BASE_DNA_DIR" -type f -path "*dna_*/*_sorted.bam" | sort)
mapfile -t RNA_FILES < <(find "$BASE_RNA_DIR" -type f -path "*rna_*/*_sorted.bam" | sort)


total_dna=${#DNA_FILES[@]}
total_rna=${#RNA_FILES[@]}

echo "Found $total_dna DNA BAMs and $total_rna RNA BAMs."
echo "Detected $total_pairs valid matched pairs."

if [[ $total_dna -eq 0 && $total_rna -eq 0 ]]; then
    echo "No BAM files found!"
    exit 1
fi

# =========================================================
# Select pairs
# =========================================================

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

# Проверяем и создаём RAW_CONTACTS_DIR, если нужно
if [[ -n "$RAW_CONTACTS_DIR" ]]; then
    mkdir -p "$RAW_CONTACTS_DIR"
else
    echo "Error: --raw-contacts-dir must be specified."
    exit 1
fi
    
    
# =========================================================
# Обработка DNA-файлов (тип эксперимента OTA_SE)
# =========================================================
if [[ $total_dna -gt 0 ]]; then
    echo "Processing DNA files (OTA_SE)...${CONTACT_SCRIPT}"
    for dna_file in "${DNA_FILES[@]}"; do
        # Формируем уникальное имя выходной подпапки на основе имени директории с BAM
        dna_dir=$(basename "$(dirname "$dna_file")")
        out_path="$RAW_CONTACTS_DIR/$dna_dir"
        mkdir -p "$out_path"
        cd "$out_path" || exit

        echo "  DNA: $dna_file -> $out_path"
        python "$CONTACT_SCRIPT" \
            -r1 "$dna_file" \
            -mr "$DNA_MODE" \
            -e "OTA_SE" \
            -p "raw_contacts"
    done
fi

# =========================================================
# Обработка RNA-файлов (тип эксперимента RNA_SEQ_SE)
# =========================================================
if [[ $total_rna -gt 0 ]]; then
    echo "Processing RNA files (RNA_SEQ_SE)..."
    for rna_file in "${RNA_FILES[@]}"; do
        rna_dir=$(basename "$(dirname "$rna_file")")
        out_path="$RAW_CONTACTS_DIR/$rna_dir"
        mkdir -p "$out_path"
        cd "$out_path" || exit

        echo "  RNA: $rna_file -> $out_path"
        python "$CONTACT_SCRIPT" \
            -r1 "$rna_file" \
            -mr "$RNA_MODE" \
            -e "RNA_SEQ_SE" \
            -p "raw_contacts"
    done
fi
    
echo "All processing complete."

