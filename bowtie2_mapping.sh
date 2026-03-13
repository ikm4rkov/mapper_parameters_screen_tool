set -euo pipefail

usage() {
  echo "Usage: $0 -f <fastq_list> -r <star_index> -t <threads> -o <output_path> -b <star_path>"
  echo ""
  echo "Required:"
  echo "  -f   Comma-separated FASTQ list (absolute paths)"
  echo "  -r   STAR index basename"
  echo "  -t   Threads"
  echo "  -o   Output directory"
  echo "  -b   STAR installation path"
  echo ""
  exit 1
}

# ==========================================
# Parameter ranges  (UNCHANGED)
# Doing ~500 values to intersect about couple of days (25*25 was 6 hours, now having ~100-200 faster it is about 250-500 params same time)
# ==========================================
D_VALUES="5 15 20 " # 10 25
R_VALUES="1 3 4" # 2
N_VALUES="0 1"
L_VALUES="10 18 25" # 20 24
I_VALUES="S,1,1.15 S,1,0.50" # S,1,0.30 S,0.5,1.50, S,0,2.50
# MA_VALUES="1 2 3" only matter in --local mode
MP_VALUES="4,2 6,2" #"6,3 8,2"
SCORE_MIN="L,-1.0,-0.6 L,3.0,-0.6 L,2.0,-0.8" # L,0,-0.6 L,2,-0.6 L,1.0,-0.8 L,1,-0.6, L,4.0,-0.5

# ==========================================
# Argument parsing
# ==========================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) INPUT_FASTQS="$2"; shift 2 ;;
    -r) REFERENCE="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -o) OUTPUT_PATH="$2"; shift 2 ;;
    -b) BOWTIE2="$2"; shift 2 ;;
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
[[ -z "$BOWTIE2" ]]  && { echo "Error: -b is required"; usage; }

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

   for D in $D_VALUES; do
     for R in $R_VALUES; do
       for N in $N_VALUES; do
         for L in $L_VALUES; do
          for i in $I_VALUES; do
            # for ma in $MA_VALUES; do
              for mp in $MP_VALUES; do
                 for score_min in $SCORE_MIN; do

                  # Build params
                  BOWTIE2_PARAMS="-D $D -R $R -N $N -L $L -i $i --mp $mp --score-min $score_min"
                  #BOWTIE2_PARAMS="--mp $mp"
                  BOWTIE2_PARAMS="--end-to-end -k 10 $BOWTIE2_PARAMS"

                  # Create directory with counter
                  MAKE_DIR="${OUTPUT_PATH}/${PART}_bowtie2_${i}"
                  mkdir -p "$MAKE_DIR"
                  echo -e "${PART}_bowtie2_${i}\t$BOWTIE2_PARAMS" >> "$WORK_DIR/parameters_bowtie2.log"
                  i=$((i + 1))
                  # PARAMS_SUFFIX=$(echo "$BOWTIE2_PARAMS" | sed 's/-//g' | sed 's/ /_/g')
                  # MAKE_DIR="${OUTPUT_PATH}/${PART}_output_${PARAMS_SUFFIX}"
                  # mkdir -p "$MAKE_DIR"

                  OUTPUT_BAM="${MAKE_DIR}/output_sorted.bam"
                  LOG_FILE="${MAKE_DIR}/alignment.log"

                  echo "Running BOWTIE2 with params: $BOWTIE2_PARAMS" 
                  $BOWTIE2/bowtie2 $BOWTIE2_PARAMS -x "$REFERENCE" -U "$INPUT_FASTQ" -p "$THREADS" 2>"$LOG_FILE" | samtools view -hb - 2>>"$LOG_FILE" | samtools sort -n -o "$OUTPUT_BAM" 2>>"$LOG_FILE"

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
            # done
          done
         done
       done
     done
   done
done

echo "Script completed successfully."

