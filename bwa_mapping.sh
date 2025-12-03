#!/bin/bash
# BWA-MEM Parameter Sweep Script
# ================================



# Аргументы: -f <fastq_list> -r <reference> -t <threads> -o <output_path> [-m <1|0>] -b <bwa_mem_path>
while getopts "f:r:t:o:m:b:" opt; do
  case "$opt" in
    f) INPUT_FASTQS="$OPTARG" ;;
    r) REFERENCE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    o) OUTPUT_PATH="$OPTARG" ;;
    m) MARGI_MODE="$OPTARG" ;;
    b) BWA_MEM="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done


if [ -z "$INPUT_FASTQS" ] || [ -z "$REFERENCE" ] || [ -z "$THREADS" ] || [ -z "$OUTPUT_PATH" ] || [ -z "$BWA_MEM" ]; then
  echo "Usage: $0 -f <fastq_list> -r <reference> -t <threads> -o <output_path> [-m <1|0>] -b <bwa_mem_path>"
  exit 1
fi

# ================================
# Parameter ranges
# ================================
T_VALUES="1" #10 20 30
K_VALUES="15 19" #5 10 
L_VALUES="4" #"3 4 5"
# B_VALUES="2 4"
# O_VALUES="4 6"
#W_VALUES="1 25 50 100"


# Make output dir (absolute)
mkdir -p "$OUTPUT_PATH"

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

  # sanity check FASTQ file exists
  if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Warning: FASTQ file not found: "$INPUT_FASTQ" - skipping" >&2
    continue
  fi

  PART=$(basename "$INPUT_FASTQ" | cut -d'.' -f2)
  
  for T in $T_VALUES; do
    for K in $K_VALUES; do
      for L in $L_VALUES; do
       # for W in $W_VALUES; do

          BWA_PARAMS="-Y -k $K -T $T -L $L"
          if [ "$MARGI_MODE" -eq 1 ]; then
            BWA_PARAMS="$BWA_PARAMS -5"
          fi

          PARAMS_SUFFIX=$(echo "$BWA_PARAMS" | sed 's/-//g' | sed 's/ /_/g')
          MAKE_DIR="${OUTPUT_PATH}/${PART}_output_${PARAMS_SUFFIX}"
          mkdir -p "$MAKE_DIR"

          OUTPUT_BAM="${MAKE_DIR}/output.bam"
          LOG_FILE="${MAKE_DIR}/alignment.log"

          echo "Running BWA with params: $BWA_PARAMS" | tee -a "$LOG_FILE"

          if ! $BWA_MEM/bwa mem -t "$THREADS" $BWA_PARAMS "$REFERENCE" "$INPUT_FASTQ" \
               2>"$LOG_FILE" | samtools view -hb - 2>>"$LOG_FILE" | samtools sort -n -o "$OUTPUT_BAM" 2>>"$LOG_FILE"; then
            echo "Error: mapping/sorting failed for $INPUT_FASTQ with params: $BWA_PARAMS" | tee -a "$LOG_FILE" >&2
            continue
          fi

          # verify bam was created
          if [ ! -f "$OUTPUT_BAM" ]; then
            echo "Error: output BAM not created: $OUTPUT_BAM" | tee -a "$LOG_FILE" >&2
            continue
          fi

        #done
      done
    done
  done
done

echo "Script completed successfully."

