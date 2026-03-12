#!/usr/bin/env bash
set -euo pipefail

###############################################
# Examples:
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode BWA --dna-mode BWA
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode STAR --dna-mode STAR
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode HISAT2 --dna-mode HISAT2
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode BOWTIE2 --dna-mode BOWTIE2
#
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode STAR --dna-mode BWA
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode HISAT2 --dna-mode BWA
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode HISAT2 --dna-mode BOWTIE2
# ./master_until_cigar.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --rna-mode STAR --dna-mode BOWTIE2
###############################################

###############################################
# Usage function
###############################################
# usage() {
#     echo "Usage: $0 <config_file> -b <base_dir> --mode <STAR|HISAT2|BWA> [options]"
#     echo
#     echo "--mode                        Mapper: STAR, HISAT2, BWA (подходит и для BOWTIE2)"
#     echo
#     echo "Options:"
#     # echo "  -b, --base-dir DIR          Base directory for BAM/SAM processing (required)"
#     echo "  --sample-contacts 0|1       Whether to run sample contacts script (default 0)"
#     echo "  --contacts-mode all|one              Mode for contacts generation script (default all)"
#     echo "  --pairs FILE                Pairs TSV file for contacts script (default pairs.tsv)"
#     echo "  --mapper-branch native|collab  Branch for script #2 (default collab)"
#     echo
#     echo "Options for script #2 random mode:"
#     echo "  -c, --count NUMBER"
#     echo "  -s, --seed NUMBER"
#     echo "  -d, --different 0|1"
#     exit 1
# }

# Получаем путь до конфига как первый аргумент
if [ $# -lt 1 ]; then
    usage
fi
config_file="$1"
shift

# Check config file existence and source it (required)
if [ ! -f "$config_file" ]; then
  echo "Error: Config file '$config_file' not found. This parameter is required."
  exit 1
fi
# shellcheck source=/dev/null
source "$config_file"

###############################################
# Parse command-line arguments
###############################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        --rna-mode)             RNA_MODE="$2"; shift 2 ;;
        --dna-mode)             DNA_MODE="$2"; shift 2 ;;
        # -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done


if [[ -z "$RNA_MODE" ]]; then
    echo "Error: --RNA_MODE is required."
    usage
fi

if [[ -z "$DNA_MODE" ]]; then
    echo "Error: --DNA_MODE is required."
    usage
fi

# RAW_CONTACTS_DIR="$WORK_DIR/raw_contacts/${MODE,,}"
RAW_CONTACTS_DIR="$WORK_DIR/raw_contacts"
# BASE_DIR="$WORK_DIR/mapping/${MODE,,}"
BASE_RNA_DIR="$WORK_DIR/mapping/${RNA_MODE,,}"
BASE_DNA_DIR="$WORK_DIR/mapping/${DNA_MODE,,}"

LOG_FILE="wrapper_rna_${RNA_MODE,,}__dna_${DNA_MODE,,}.log"
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

LOG_FILE="$RAW_CONTACTS_DIR/$LOG_FILE"
echo "Logging to $LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

###############################################
# Step 2: Generate contacts
###############################################
SCRIPT2_ARGS=(--base-rna-dir "$BASE_RNA_DIR" --base-dna-dir "$BASE_DNA_DIR" --raw-contacts-dir "$RAW_CONTACTS_DIR" -k "$RNA_MODE" -kk "$DNA_MODE" --contact_script "$CONTACT_SCRIPT")
./make_contacts_universal_notsample.sh "${SCRIPT2_ARGS[@]}"

FILTER_OUTPUT_DIR="$WORK_DIR/filtered_contacts"
mkdir -p "$FILTER_OUTPUT_DIR"
./calculate_CIAGAR_filter.sh -i "$RAW_CONTACTS_DIR" -o "$FILTER_OUTPUT_DIR" -s "$EDIT_DISTANCE_CIGAR_FILTER_SCRIPT" -k "$RNA_MODE" -kk "$DNA_MODE"