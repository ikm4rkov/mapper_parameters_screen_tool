#!/usr/bin/env bash

set -euo pipefail

##############################################
# Helper: print usage
##############################################
usage() {
    cat <<EOF
Usage: $0 --imargi <true|false> --mode <bowtie2|star|hisat2|bwa> [OPTIONS]

Common options:
  --imargi        true or false
  --mode          Mapper: bowtie2, star, hisat2, bwa

Bowtie2 / STAR options:
  --input         Input FASTQ file
  --part          dna or rna
  --seed          Random seed
  --config-log    Optional log name
  --index         Path to genome index
  --threads       Number of threads

BWA / HISAT2 options:
  -f              FASTQ inputs (comma-separated)
  -r              Reference genome index
  -t              Number of threads
  -o              Output directory
EOF
    exit 1
}

##############################################
# Argument parsing
##############################################

# Defaults
IMARGI=""
MODE=""
INPUT=""
PART=""
SEED=""
CONFIG_LOG=""
INDEX=""
THREADS=""
F_ARG=""
R_ARG=""
T_ARG=""
O_ARG=""

# Manual argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        --imargi)      IMARGI="$2"; shift 2 ;;
        --mode)        MODE="$2"; shift 2 ;;
        --input)       INPUT="$2"; shift 2 ;;
        --part)        PART="$2"; shift 2 ;;
        --seed)        SEED="$2"; shift 2 ;;
        --config-log)  CONFIG_LOG="$2"; shift 2 ;;
        --index)       INDEX="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        -f)            F_ARG="$2"; shift 2 ;;
        -r)            R_ARG="$2"; shift 2 ;;
        -t)            T_ARG="$2"; shift 2 ;;
        -o)            O_ARG="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *) echo "Unknown argument: $1" ; usage ;;
    esac
done

##############################################
# Validate --imargi
##############################################
if [[ "$IMARGI" == "true" ]]; then
    echo "not implemented yet"
    exit 0
elif [[ "$IMARGI" != "false" ]]; then
    echo "Error: --imargi must be true or false."
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
        CMD=( python3 "$WORKER_SCRIPT"
              --input "$INPUT"
              --part "$PART"
              --seed "$SEED"
              --index "$INDEX"
              --threads "$THREADS" )
        [[ -n "$CONFIG_LOG" ]] && CMD+=( --config-log "$CONFIG_LOG" )
        ;;

    star)
        WORKER_SCRIPT="${SCRIPT_DIR}/star_montecarlo.py"
        CMD=( python3 "$WORKER_SCRIPT"
              --input "$INPUT"
              --part "$PART"
              --seed "$SEED"
              --index "$INDEX"
              --threads "$THREADS" )
        [[ -n "$CONFIG_LOG" ]] && CMD+=( --config-log "$CONFIG_LOG" )
        ;;

    bwa)
        WORKER_SCRIPT="${SCRIPT_DIR}/bwa_mapping.sh"
        CMD=( bash "$WORKER_SCRIPT"
              -f "$F_ARG"
              -r "$R_ARG"
              -t "$T_ARG"
              -o "$O_ARG" )
        ;;

    hisat2)
        WORKER_SCRIPT="${SCRIPT_DIR}/hisat2_mapping.sh"
        CMD=( bash "$WORKER_SCRIPT"
              -f "$F_ARG"
              -r "$R_ARG"
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

