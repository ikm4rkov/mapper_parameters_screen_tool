et -euo pipefail

##############################################
# Usage message
##############################################
usage() {
        cat <<EOF
Usage: $0 <config_file> --mode <bowtie2|star|hisat2|bwa> [OPTIONS]

config_file: path to mappers.conf

Common parameters:
  --mode              Mapper to run
  --imargi            1 or 0
  --part              dna/rna (bowtie2 + star)
  --seed              Random seed (bowtie2 + star)
  --config-log        Name for config log (bowtie2 + star)
  --threads           Threads (bowtie2 + star)
  -t                  Threads (bwa/hisat2)

FASTQ_INPUTS and OUTPUT_PREFIX must be supplied via config file.

Reference index convention:
  \$WORK_DIR/genomes/<mapper>/\$OUTPUT_PREFIX

Output directory convention:
  \$WORK_DIR/mapping/<mapper>/
EOF
        exit 1
}

##############################################
# Parse input
##############################################
if [ $# -lt 1 ]; then usage; fi

config_file="$1"
shift

if [[ ! -f "$config_file" ]]; then
    echo "Error: config file '$config_file' not found."
    exit 1
fi

source "$config_file"

# ---- defaults ----
IMARGI="0"
MODE=""
PART=""
SEED=""
CONFIG_LOG=""
THREADS=""
T_ARG=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

##############################################
# Minimal manual parser
##############################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        --imargi)      IMARGI="$2"; shift 2 ;;
        --mode)        MODE="$2"; shift 2 ;;
        --part)        PART="$2"; shift 2 ;;
        --seed)        SEED="$2"; shift 2 ;;
        --config-log)  CONFIG_LOG="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        -t)            T_ARG="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

##############################################
# Validate FASTQ_INPUTS / OUTPUT_PREFIX
##############################################
if [[ -z "${FASTQ_INPUTS:-}" ]]; then
    echo "ERROR: FASTQ_INPUTS must be exported in the config file."
    exit 1
fi

if [[ -z "${OUTPUT_PREFIX:-}" ]]; then
    echo "ERROR: OUTPUT_PREFIX must be exported in the config file."
    exit 1
fi

##############################################
# Standardized reference and output dirs
##############################################
REFDIR="$WORK_DIR/genomes"
OUTDIR="$WORK_DIR/mapping"

##############################################
# Mode dispatch
##############################################
CMD=()
WORKER_SCRIPT=""

case "$MODE" in

################################################
# BOWTIE2 (Python randomized sampling)
################################################
    bowtie2)
        WORKER_SCRIPT="${SCRIPT_DIR}/bowtie2_montecarlo.py"

        CMD=(
            python3 "$WORKER_SCRIPT"
            --fastq-inputs "$FASTQ_INPUTS"
            --part "$PART"
            --seed "$SEED"
            --index "$REFDIR/bowtie2/$OUTPUT_PREFIX"
            --threads "$THREADS"
            --bowtie2-bin "$BOWTIE2"
        )

        [[ -n "$CONFIG_LOG" ]] && CMD+=( --config-log "$CONFIG_LOG" )
        ;;

################################################
# STAR (Python randomized sampling)
################################################
    star)
        WORKER_SCRIPT="${SCRIPT_DIR}/star_montecarlo.py"

        CMD=(
            python3 "$WORKER_SCRIPT"
            --fastq-inputs "$FASTQ_INPUTS"
            --part "$PART"
            --seed "$SEED"
            --index "$REFDIR/star/$OUTPUT_PREFIX"
            --threads "$THREADS"
            --star-bin "$STAR"
        )

        [[ -n "$CONFIG_LOG" ]] && CMD+=( --config-log "$CONFIG_LOG" )
        ;;

################################################
# BWA (shell script)
################################################
    bwa)
        WORKER_SCRIPT="${SCRIPT_DIR}/bwa_mapping.sh"

        CMD=(
            bash "$WORKER_SCRIPT"
            -f "$FASTQ_INPUTS"
            -r "$REFDIR/bwa/$OUTPUT_PREFIX"
            -t "$T_ARG"
            -o "$OUTDIR/bwa"
            -m "$IMARGI"
            -b "$BWA_MEM"
        )
        ;;

################################################
# HISAT2 (shell script)
################################################
    hisat2)
        WORKER_SCRIPT="${SCRIPT_DIR}/hisat2_mapping.sh"

        CMD=(
            bash "$WORKER_SCRIPT"
            -f "$FASTQ_INPUTS"
            -r "$REFDIR/hisat2/$OUTPUT_PREFIX"
            -t "$T_ARG"
            -o "$OUTDIR/hisat2"
            -b "$HISAT2"
        )
        ;;

    *)
        echo "Unknown mode: $MODE"
        exit 1
        ;;
esac

##############################################
# Check worker exists
##############################################
if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "ERROR: Worker script not found: $WORKER_SCRIPT"
    exit 1
fi

##############################################
# Run
##############################################
echo "Executing: ${CMD[*]}"
exec "${CMD[@]}"

