#!/bin/bash

# ==========================================
# Default configuration
# ==========================================
HISAT2="hisat2"
THREADS=15
REFERENCE="/home/ryabykh2018/iMARGI/alignment_tools/data/Sus_scrofa/genome/hisat2/Sus_scrofa.Sscrofa11.1_index"

INPUT_FASTQ_PATH="/home/ryabykh2018/iMARGI/alignment_tools/data/Sus_scrofa/subsample_1M"
OUTPUT_PATH="/home/ryabykh2018/iMARGI/alignment_tools/data/Sus_scrofa/subsample_1M/hisat2_test"

MARGI_MODE=0   # default: off

# ==========================================
# Parameter ranges
# ==========================================
IGNORE_QUALS_VALUES="--ignore-quals ''"
NO_TEMPLATELEN_ADJUSMENT_VALUES="--no-templatelen-adjustment ''"
NON_DETERMINISTIC_VALUES="--non-deterministic ''"
PEN_NONCANSPLICE_VALUES="12"
MAX_SEEDS_VALUES="5 200"
SCORE_MIN_VALUES="L,0,-0.4"
MP_VALUES="6,2"
RDG_VALUES="5,3"
RFG_VALUES="5,3"

# ==========================================
# Usage function
# ==========================================
usage() {
  echo "Usage: $0 -f <fastq_list> [-r <reference>] [-t <threads>] [-o <output_path>] [-m <1|0>]"
  echo ""
  echo "  -f   Comma-separated FASTQ list"
  echo "  -r   HISAT2 index basename"
  echo "  -t   Threads (default 15)"
  echo "  -o   Output directory"
  echo "  -m   MARGI mode: 1 adds -SP5M, 0 disables (default 0)"
  echo "       Long option: --margi 0|1"
  exit 1
}

# ==========================================
# Argument parsing
# ==========================================
# Support long option --margi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) INPUT_FASTQS="$2"; shift 2 ;;
    -r) REFERENCE="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -o) OUTPUT_PATH="$2"; shift 2 ;;
    -m|--margi)
        MARGI_MODE="$2"
        shift 2
        ;;
    -*)
        echo "Unknown option: $1"
        usage
        ;;
    *)
        break
        ;;
  esac
done

if [ -z "$INPUT_FASTQS" ]; then
  echo "Error: please provide FASTQ list with -f"
  usage
fi

# Validate MARGI_MODE
if [[ "$MARGI_MODE" != "0" && "$MARGI_MODE" != "1" ]]; then
  echo "Error: -m / --margi must be 0 or 1"
  exit 1
fi

if [[ "$MARGI_MODE" -eq 1 ]]; then
  echo "MARGI mode ON → adding -SP5M to HISAT2 command"
else
  echo "MARGI mode OFF"
fi

# ==========================================
# Preparations
# ==========================================
mkdir -p "$OUTPUT_PATH"

if ! command -v "$HISAT2" >/dev/null 2>&1; then
  echo "Error: HISAT2 not found: $HISAT2" >&2
  exit 1
fi

trim() {
  echo "$1" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

# ==========================================
# Main loop
# ==========================================
IFS=',' read -ra FASTQ_ARRAY <<< "$INPUT_FASTQS"
for raw in "${FASTQ_ARRAY[@]}"; do
  INPUT_FASTQ=$(trim "$raw")

  if [[ "$INPUT_FASTQ" = /* ]]; then
    FASTQ_FULLPATH="$INPUT_FASTQ"
  else
    FASTQ_FULLPATH="${INPUT_FASTQ_PATH%/}/$INPUT_FASTQ"
  fi

  if [ ! -f "$FASTQ_FULLPATH" ]; then
    echo "Warning: FASTQ not found: $FASTQ_FULLPATH — skipping" >&2
    continue
  fi

  PART=$(basename "$INPUT_FASTQ" | cut -d'.' -f2)
  N_READS_LOG="${OUTPUT_PATH}/n_reads_${PART}.log"
  echo -e "params\tN_uniq_mapped_reads\tTotal_unique_read_IDs" > "$N_READS_LOG"

  for ignore_quals in $IGNORE_QUALS_VALUES; do
    for no_templatelen in $NO_TEMPLATELEN_ADJUSMENT_VALUES; do
      for nondeterministic in $NON_DETERMINISTIC_VALUES; do
        for pen_noncansplice in $PEN_NONCANSPLICE_VALUES; do
          for max_seeds in $MAX_SEEDS_VALUES; do
            for score_min in $SCORE_MIN_VALUES; do
              for mp in $MP_VALUES; do
                for rdg in $RDG_VALUES; do
                  for rfg in $RFG_VALUES; do

                    HISAT2_PARAMS="$(echo $ignore_quals | tr -d "'") $(echo $no_templatelen | tr -d "'") $(echo $nondeterministic | tr -d "'") --pen-noncansplice $pen_noncansplice --max-seeds $max_seeds --score-min $score_min --mp $mp --rdg $rdg --rfg $rfg"

                    # Add MARGI parameter
                    if [[ "$MARGI_MODE" -eq 1 ]]; then
                      HISAT2_PARAMS="$HISAT2_PARAMS -SP5M"
                    fi

                    PARAMS_SUFFIX=$(echo "$HISAT2_PARAMS" | sed 's/-//g' | sed 's/ /_/g')
                    MAKE_DIR="${OUTPUT_PATH}/${PART}_output_${PARAMS_SUFFIX}"
                    mkdir -p "$MAKE_DIR"

                    OUTPUT_BAM="${MAKE_DIR}/output.bam"
                    LOG_FILE="${MAKE_DIR}/hisat2_summary.txt"

                    echo "Running HISAT2 with params: $HISAT2_PARAMS" | tee -a "$LOG_FILE"

                    if ! "$HISAT2" -x "$REFERENCE" -U "$FASTQ_FULLPATH" -k 100 --dta \
                      $HISAT2_PARAMS -p "$THREADS" --summary-file "$LOG_FILE" 2>>"$LOG_FILE" \
                      | samtools view -@ "$THREADS" -Sb - 2>>"$LOG_FILE" \
                      | samtools sort -n -o "$OUTPUT_BAM" 2>>"$LOG_FILE"; then

                      echo "Error: HISAT2 failed for $FASTQ_FULLPATH" | tee -a "$LOG_FILE" >&2
                      continue
                    fi

                    if [ ! -f "$OUTPUT_BAM" ]; then
                      echo "Error: BAM not created" | tee -a "$LOG_FILE" >&2
                      continue
                    fi

                    COUNT_UNIQ_MAPPED=$(samtools view -F 256 -F 4 -e '[NH]==1' "$OUTPUT_BAM" | wc -l)
                    COUNT_TOTAL_UNIQ_IDS=$(samtools view "$OUTPUT_BAM" | cut -f1 | sort | uniq | wc -l)

                    echo -e "${PARAMS_SUFFIX}\t${COUNT_UNIQ_MAPPED}\t${COUNT_TOTAL_UNIQ_IDS}" >> "$N_READS_LOG"

                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

echo "Script completed successfully."
   
