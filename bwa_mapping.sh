#!/bin/bash: indent: command not found-MEM2 Parameter Sweep Script
# ================================

BWA_MEM2="bwa-mem2"
THREADS=8
REFERENCE="/home/ryabykh2018/iMARGI/alignment_tools/data/Sus_scrofa/genome/bwa_mem2/Sus_scrofa.Sscrofa11.1.dna.toplevel"
OUTPUT_PATH="/home/ryabykh2018/iMARGI/alignment_tools/data/Sus_scrofa/subsample_1M/bwa_mem2_test"
INPUT_FASTQ_PATH="/home/ryabykh2018/iMARGI/alignment_tools/data/Sus_scrofa/subsample_1M"
MARGI_MODE=0  # Default: off

# ================================
# Parameter ranges
# ================================
T_VALUES="1 10 20 30"
K_VALUES="5 10 15 19"
L_VALUES="3 4 5"
# B_VALUES="2 4"
# O_VALUES="4 6"
#W_VALUES="1 25 50 100"

# ================================
# Usage
# ================================
usage() {
  echo "Usage: $0 -f <fastq_list> [-r <reference>] [-t <threads>] [-o <output_path>] [-m <1|0>]"
  echo "  -f fastq_list : Comma-separated list of FASTQ files (absolute paths or basenames)"
  echo "  -r reference  : Reference genome basename or path (optional)"
  echo "  -t threads    : Threads for BWA (default 8)"
  echo "  -o output     : Output directory (optional)"
  echo "  -m MARGI mode : 1 to enable '-SP5M', 0 to disable (default 0)"
  exit 1
}

# ================================
# Parse args
# ================================
while getopts "f:r:t:o:m:" opt; do
  case "$opt" in
    f) INPUT_FASTQS="$OPTARG" ;;
    r) REFERENCE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    o) OUTPUT_PATH="$OPTARG" ;;
    m) MARGI_MODE="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

if [ -z "$INPUT_FASTQS" ]; then
  echo "Error: please provide FASTQ list using -f"
  usage
fi

# Make output dir (absolute)
mkdir -p "$OUTPUT_PATH"

# Optional: check bwa executable
if ! command -v "$BWA_MEM2" >/dev/null 2>&1; then
  echo "Error: bwa executable '$BWA_MEM2' not found in PATH. Set BWA_MEM2 to full path." >&2
  exit 1
fi

# ================================
# Helper: trim whitespace (leading/trailing)
# ================================
trim() {
  # usage: trimmed=$(trim "$var")
  echo "$1" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

# ================================
# Main loop: split comma-separated list safely
# ================================
IFS=',' read -ra FASTQ_ARRAY <<< "$INPUT_FASTQS"

for raw in "${FASTQ_ARRAY[@]}"; do
  INPUT_FASTQ=$(trim "$raw")

  # If the provided FASTQ is an absolute path (starts with '/'), use it as-is.
  # Otherwise prepend INPUT_FASTQ_PATH.
  if [[ "$INPUT_FASTQ" = /* ]]; then
    FASTQ_FULLPATH="$INPUT_FASTQ"
  else
    FASTQ_FULLPATH="${INPUT_FASTQ_PATH%/}/$INPUT_FASTQ"
  fi

  # sanity check FASTQ file exists
  if [ ! -f "$FASTQ_FULLPATH" ]; then
    echo "Warning: FASTQ file not found: $FASTQ_FULLPATH - skipping" >&2
    continue
  fi

  PART=$(basename "$INPUT_FASTQ" | cut -d'.' -f2)
  N_READS_LOG="${OUTPUT_PATH}/n_reads_${PART}.log"
  echo -e "params\tN_uniq_mapped_reads\tTotal_unique_read_IDs" > "$N_READS_LOG"

  for T in $T_VALUES; do
    for K in $K_VALUES; do
      for L in $L_VALUES; do
       # for W in $W_VALUES; do

          BWA_PARAMS="-k $K -T $T -L $L"
          if [ "$MARGI_MODE" -eq 1 ]; then
            BWA_PARAMS="$BWA_PARAMS -SP5M"
          fi

          PARAMS_SUFFIX=$(echo "$BWA_PARAMS" | sed 's/-//g' | sed 's/ /_/g')
          MAKE_DIR="${OUTPUT_PATH}/${PART}_output_${PARAMS_SUFFIX}"
          mkdir -p "$MAKE_DIR"

          OUTPUT_BAM="${MAKE_DIR}/output.bam"
          LOG_FILE="${MAKE_DIR}/alignment.log"

          echo "Running BWA with params: $BWA_PARAMS" | tee -a "$LOG_FILE"

          # Run BWA -> samtools pipeline with error checking.
          if ! $BWA_MEM2 mem -t "$THREADS" $BWA_PARAMS "$REFERENCE" "$FASTQ_FULLPATH" \
               2>>"$LOG_FILE" | samtools view -hb - 2>>"$LOG_FILE" | samtools sort -n -o "$OUTPUT_BAM" 2>>"$LOG_FILE"; then
            echo "Error: mapping/sorting failed for $FASTQ_FULLPATH with params: $BWA_PARAMS" | tee -a "$LOG_FILE" >&2
            continue
          fi

          # verify bam was created
          if [ ! -f "$OUTPUT_BAM" ]; then
            echo "Error: output BAM not created: $OUTPUT_BAM" | tee -a "$LOG_FILE" >&2
            continue
          fi

          COUNT_UNIQ_MAPPED=$(samtools view -F 4 "$OUTPUT_BAM" | awk '!/SA:Z:/ && !/XA:Z:/' | wc -l)
          COUNT_TOTAL_UNIQ_IDS=$(samtools view "$OUTPUT_BAM" | cut -f1 | sort | uniq | wc -l)

          echo -e "${PARAMS_SUFFIX}\t${COUNT_UNIQ_MAPPED}\t${COUNT_TOTAL_UNIQ_IDS}" >> "$N_READS_LOG"
        #done
      done
    done
  done
done

echo "Script completed successfully."

