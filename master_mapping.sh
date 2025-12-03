#!/usr/bin/env bash

set -euo pipefail

##############################################
# Helper: print usage
##############################################
usage() {
        cat <<EOF
Usage: $0 <config_file> --imargi <1|0> --mode <bowtie2|star|hisat2|bwa> [OPTIONS]

config_file: path to the config file (mappers.conf), required

Common options:
    --imargi        1 (true) or 0 (false)
    --mode          Mapper: bowtie2, star, hisat2, bwa

Bowtie2 / STAR options:
    --part          dna or rna
    --seed          Random seed
    --config-log    Optional log name
    --threads       Number of threads

BWA / HISAT2 options:
    -t              Number of threads
    # -o              Output directory  # OUTPUT_PATH будет задан автоматически для bwa

FASTQ_INPUTS импортируется из конфига для всех режимов и не задается через параметры.
Для bwa mode, если FASTQ_INPUTS и OUTPUT_PREFIX заданы в config, они будут использованы автоматически.
Reference genome index будет, например, "$WORK_DIR/genomes/bwa/$OUTPUT_PREFIX".
EOF
        exit 1
}

##############################################
# Examples:
# ./master_mapping.sh /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf --mode bwa -t 4
##############################################

##############################################
# Argument parsing
##############################################

# Defaults
IMARGI="0"
MODE=""
PART=""
SEED=""
CONFIG_LOG=""
THREADS=""
T_ARG=""
O_ARG=""


# Получаем путь до конфига как первый аргумент
if [ $# -lt 1 ]; then
    usage
fi
config_file="$1"
shift

if [ ! -f "$config_file" ]; then
    echo "Error: config file '$config_file' not found."
    exit 1
fi
source "$config_file"

# Manual argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --imargi)      IMARGI="$2"; shift 2 ;;
        --mode)        MODE="$2"; shift 2 ;;
        --part)        PART="$2"; shift 2 ;;
        --seed)        SEED="$2"; shift 2 ;;
        --config-log)  CONFIG_LOG="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        -t)            T_ARG="$2"; shift 2 ;;
        -o)            O_ARG="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *) echo "Unknown argument: $1" ; usage ;;
    esac
done

##############################################
# Validate --imargi
##############################################
if [[ "$IMARGI" == "1" ]]; then
    echo "not implemented yet"
    exit 0
elif [[ "$IMARGI" != "0" ]]; then
    echo "Error: --imargi must be 1 or 0."
    exit 1
fi

##############################################
# Determine directory of this script
##############################################
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

##############################################
# Mode dispatch
##############################################
CMD=()
WORKER_SCRIPT=""

case "$MODE" in
    bowtie2)
        WORKER_SCRIPT="${SCRIPT_DIR}/bowtie2_montecarlo.py"
        # FASTQ_INPUTS импортируется из конфига
        CMD=( python3 "$WORKER_SCRIPT"
              --fastq-inputs "$FASTQ_INPUTS"
              --part "$PART"
              --seed "$SEED"
              --index "$WORK_DIR/genomes/bwa/$OUTPUT_PREFIX"
              --threads "$THREADS" )
        [[ -n "$CONFIG_LOG" ]] && CMD+=( --config-log "$CONFIG_LOG" )
        ;;

    star)
        WORKER_SCRIPT="${SCRIPT_DIR}/star_montecarlo.py"
        CMD=( python3 "$WORKER_SCRIPT"
              --fastq-inputs "$FASTQ_INPUTS"
              --part "$PART"
              --seed "$SEED"
              --index "$WORK_DIR/genomes/star"
              --threads "$THREADS" )
        [[ -n "$CONFIG_LOG" ]] && CMD+=( --config-log "$CONFIG_LOG" )
        ;;

    bwa)
        WORKER_SCRIPT="${SCRIPT_DIR}/bwa_mapping.sh"
        CMD=( bash "$WORKER_SCRIPT"
              -f "$FASTQ_INPUTS"
              -r "$WORK_DIR/genomes/bwa/$OUTPUT_PREFIX"
              -t "$T_ARG"
              -o "$WORK_DIR/mapping/bwa"
              -m "$IMARGI"
              -b "$BWA_MEM" )
        ;;

    hisat2)
        WORKER_SCRIPT="${SCRIPT_DIR}/hisat2_mapping.sh"
        CMD=( bash "$WORKER_SCRIPT"
              -f "$FASTQ_INPUTS"
              -r "$WORK_DIR/genomes/hisat2/$OUTPUT_PREFIX"
              -t "$T_ARG"
              -o "$O_ARG" )
        ;;

    *)
        echo "Unknown mode: $MODE"
        exit 1
        ;;
esac

##############################################
# Check that worker script exists
##############################################
if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "Error: Worker script not found: $WORKER_SCRIPT"
    exit 1
fi

##############################################
# Execute command
##############################################
echo "Executing: ${CMD[*]}"
exec "${CMD[@]}"


~

