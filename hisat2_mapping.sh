set -euo pipefail

# ==========================================
# Usage function
# ==========================================
usage() {
  echo "Usage: $0 -f <fastq_list> -r <hisat2_index> -t <threads> -o <output_path> -b <hisat2_path>"
  echo ""
  echo "Required:"
  echo "  -f   Comma-separated FASTQ list (absolute paths)"
  echo "  -r   HISAT2 index basename"
  echo "  -t   Threads"
  echo "  -o   Output directory"
  echo "  -b   HISAT2 installation path"
  echo ""
  exit 1
}

# ==========================================
# Parameter ranges  (UNCHANGED)
# ==========================================
# IGNORE_QUALS_VALUES="--ignore-quals" #"--ignore-quals ''"
# NO_TEMPLATELEN_ADJUSMENT_VALUES="--no-templatelen-adjustment" #"--no-templatelen-adjustment ''"
# NON_DETERMINISTIC_VALUES="--non-deterministic" #"--non-deterministic ''"
PEN_NONCANSPLICE_VALUES="12" # "12 1000000"
MAX_SEEDS_VALUES="200" # "5 200 2000"
SCORE_MIN_VALUES="L,0,-0.4" # "L,0,-0.4 L,0,-0.2"
MP_VALUES="6,2 3,1"
RDG_VALUES="5,3 3,1"
RFG_VALUES="5,3 3,1"

# ==========================================
# Argument parsing
# ==========================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) INPUT_FASTQS="$2"; shift 2 ;;
    -r) REFERENCE="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -o) OUTPUT_PATH="$2"; shift 2 ;;
    -b) HISAT2="$2"; shift 2 ;;
    -n) N_READS_INPUT="$2"; shift 2 ;;
    -w) WORK_DIR="$2"; shift 2 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# ==========================================
# Validate required arguments
# ==========================================
[[ -z "$INPUT_FASTQS" ]] && { echo "Error: -f is required"; usage; }
[[ -z "$REFERENCE" ]]    && { echo "Error: -r is required"; usage; }
[[ -z "$THREADS" ]]      && { echo "Error: -t is required"; usage; }
[[ -z "$OUTPUT_PATH" ]]  && { echo "Error: -o is required"; usage; }
[[ -z "$HISAT2" ]]  && { echo "Error: -b is required"; usage; }

# ==========================================
# Preparations
# ==========================================
mkdir -p "$OUTPUT_PATH"

trim() {
  echo "$1" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

# ==========================================
# Main loop
# ==========================================
IFS=',' read -ra FASTQ_ARRAY <<< "$INPUT_FASTQS"

i=1

for raw in "${FASTQ_ARRAY[@]}"; do
  INPUT_FASTQ=$(trim "$raw")

  # FASTQ must exist exactly as provided
  if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Warning: FASTQ not found: $INPUT_FASTQ — skipping" >&2
    continue
  fi

  PART=$(basename "$INPUT_FASTQ" | cut -d'.' -f2)

  # for ignore_quals in $IGNORE_QUALS_VALUES; do
  #   for no_templatelen in $NO_TEMPLATELEN_ADJUSMENT_VALUES; do
  #     for nondeterministic in $NON_DETERMINISTIC_VALUES; do
         for pen_noncansplice in $PEN_NONCANSPLICE_VALUES; do
          for max_seeds in $MAX_SEEDS_VALUES; do
             for score_min in $SCORE_MIN_VALUES; do
               for mp in $MP_VALUES; do
                 for rdg in $RDG_VALUES; do
                   for rfg in $RFG_VALUES; do

                    # Build params
                    # HISAT2_PARAMS="$(echo $ignore_quals | tr -d "'") $(echo $no_templatelen | tr -d "'") $(echo $nondeterministic | tr -d "'") --pen-noncansplice $pen_noncansplice --max-seeds $max_seeds --score-min $score_min --mp $mp --rdg $rdg --rfg $rfg"
                    HISAT2_PARAMS="--pen-noncansplice $pen_noncansplice --max-seeds $max_seeds --score-min $score_min --mp $mp --rdg $rdg --rfg $rfg"

                    # Добавляем --dta-cufflinks для PART==rna, иначе --no-spliced-alignment
                    if [[ "$PART" == "rna" ]]; then
                      HISAT2_PARAMS="$HISAT2_PARAMS --dta-cufflinks"
                    else
                      HISAT2_PARAMS="$HISAT2_PARAMS --no-spliced-alignment"
                    fi

                    # Create directory with counter
                    MAKE_DIR="${OUTPUT_PATH}/${PART}_hisat2_${i}"
                    mkdir -p "$MAKE_DIR"
                    echo -e "${PART}_hisat2_${i}\t$HISAT2_PARAMS" >> "$WORK_DIR/parameters_hisat2.log"
                    i=$((i + 1))
                    # PARAMS_SUFFIX=$(echo "$HISAT2_PARAMS" | sed 's/-//g' | sed 's/ /_/g')
                    # MAKE_DIR="${OUTPUT_PATH}/${PART}_output_${PARAMS_SUFFIX}"
                    # mkdir -p "$MAKE_DIR"

                    OUTPUT_BAM="${MAKE_DIR}/output_sorted.bam"
                    LOG_FILE="${MAKE_DIR}/hisat2_summary.txt"

                    echo "Running HISAT2 with params: $HISAT2_PARAMS" | tee -a "$LOG_FILE"
                    
                    $HISAT2/hisat2 -x "$REFERENCE" -U "$INPUT_FASTQ" -k 100 $HISAT2_PARAMS -p "$THREADS" \
                      --summary-file "$LOG_FILE" | samtools view -Sb - | samtools sort -n -o "$OUTPUT_BAM"

                    if [ ! -f "$OUTPUT_BAM" ]; then
                      echo "Error: BAM not created" | tee -a "$LOG_FILE" >&2
                      continue
                    fi

                    # Check if the number of unique reads matches N_READS_INPUT
                    if [[ $(samtools view "$OUTPUT_BAM" | cut -f1 | sort | uniq | wc -l) -ne "$N_READS_INPUT" ]]; then
                      echo "Error: The number of unique reads does not match N_READS_INPUT ($N_READS_INPUT)." | tee -a "$LOG_FILE" >&2
                      continue
                    fi

                  done
                 done
               done
             done
           done
         done
  #     done
  #   done
  # done
done

echo "Script completed successfully."

