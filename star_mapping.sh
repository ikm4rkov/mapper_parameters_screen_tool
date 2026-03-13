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
# Parameter ranges  (UNCHANGED) ~500 values
# ==========================================
SCORE_GAP="-5 0" #-10 5 10
SCORE_GAP_NONCAN="-8 -4" #-16 -12 0
# SCORE_GAP_GCAG="-8 -6 -4 -2 0"
# SCORE_GAP_ATAC="-16 -12 -8 -4 0"
SCORE_DEL_OPEN="-3 -2" #-4 -1 0
SCORE_INS_OPEN="-2 -1" #-4 -3 0
# SCORE_STITCH_SJSHIFT="-4 -2 1 3 5"
# SEED_SEARCH_START_LMAX="25 37.5 50 62.5 75"
SEED_SEARCH_START_LMAX_OVER_LREAD="0.75 1.0 1.25" #0.5 1,5
# SEED_PER_READ_NMAX="500 750 1000 1250 1500"
SEED_PER_WINDOW_NMAX="25 50 75" #37.5 62.5
# SEED_SPLIT_MIN="6 9 12 15 18"
# WIN_ANCHOR_DIST_NBINS="4 6 9 11 13"
OUT_FILTER_SCORE_MIN="0 5" #"-5 -2.5 0 2.5 5"
OUT_FILTER_MATCH_N_MIN_OVER_LREAD="0.5 0.66"

# ==========================================
# Argument parsing
# ==========================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) INPUT_FASTQS="$2"; shift 2 ;;
    -r) REFERENCE="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -o) OUTPUT_PATH="$2"; shift 2 ;;
    -b) STAR="$2"; shift 2 ;;
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
[[ -z "$STAR" ]]  && { echo "Error: -b is required"; usage; }

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

  for scoreGap in $SCORE_GAP; do
     for scoreGapNoncan in $SCORE_GAP_NONCAN; do
  #     for scoreGapGCAG in $SCORE_GAP_GCAG; do
  #       for scoreGapATAC in $SCORE_GAP_ATAC; do
           for scoreDelOpen in $SCORE_DEL_OPEN; do
             for scoreInsOpen in $SCORE_INS_OPEN; do
  #             for scoreStitchSJshift in $SCORE_STITCH_SJSHIFT; do
  #               for seedSearchStartLmax in $SEED_SEARCH_START_LMAX; do
                   for seedSearchStartLmaxOverLread in $SEED_SEARCH_START_LMAX_OVER_LREAD; do
  #                   for seedPerReadNmax in $SEED_PER_READ_NMAX; do
                       for seedPerWindowNmax in $SEED_PER_WINDOW_NMAX; do
  #                       for seedSplitMin in $SEED_SPLIT_MIN; do
  #                         for winAnchorDistNbins in $WIN_ANCHOR_DIST_NBINS; do
                            for outFilterScoreMin in $OUT_FILTER_SCORE_MIN; do
			      for outFilterMatchNminOverLread in $OUT_FILTER_MATCH_N_MIN_OVER_LREAD; do
                              
                              # Build params
                              STAR_PARAMS="--scoreGap $scoreGap --scoreGapNoncan $scoreGapNoncan --scoreDelOpen $scoreDelOpen --scoreInsOpen $scoreInsOpen --seedSearchStartLmaxOverLread $seedSearchStartLmaxOverLread --seedPerWindowNmax $seedPerWindowNmax --seedSplitMin $seedSplitMin --outFilterScoreMin $outFilterScoreMin --outFilterMatchNminOverLread $outFilterMatchNminOverLread"
                              # STAR_PARAMS="--outFilterScoreMin $outFilterScoreMin"
                              STAR_PARAMS="$STAR_PARAMS --outSAMattributes NM --outSAMunmapped Within"
                              
                              if [[ "$PART" == "dna" ]]; then
                                STAR_PARAMS="$STAR_PARAMS --alignIntronMax 1"  
                              fi

                              # Create directory with counter
                              MAKE_DIR="${OUTPUT_PATH}/${PART}_star_${i}"
                              mkdir -p "$MAKE_DIR"
                              echo -e "${PART}_star_${i}\t$STAR_PARAMS" >> "$WORK_DIR/parameters_star.log"
                              i=$((i + 1))
                              
                              # PARAMS_SUFFIX=$(echo "$STAR_PARAMS" | sed 's/-//g' | sed 's/ /_/g')
                              # MAKE_DIR="${OUTPUT_PATH}/${PART}_output_${PARAMS_SUFFIX}"
                              # mkdir -p "$MAKE_DIR"

                              OUTPUT_BAM="${MAKE_DIR}/output_"
                              LOG_FILE="${MAKE_DIR}/star_summary.txt"

                              echo "Running STAR with params: $STAR_PARAMS" | tee -a "$LOG_FILE"
                              $STAR/STAR --genomeDir "$REFERENCE" --readFilesIn "$INPUT_FASTQ" $STAR_PARAMS --outFileNamePrefix "$OUTPUT_BAM" --runThreadN "$THREADS"

                              SAM_FILE="${OUTPUT_BAM}Aligned.out.sam"
                              SORTED_BAM="${MAKE_DIR}/output_sorted.bam"
                              
                              if [ ! -f "$SAM_FILE" ]; then
                                echo "Error: SAM file not created: $SAM_FILE" | tee -a "$LOG_FILE" >&2
                                continue
                              fi

                              echo "Found SAM file: $SAM_FILE, running samtools sort..." | tee -a "$LOG_FILE"
                              samtools sort -n "$SAM_FILE" -o "$SORTED_BAM"
                              
                              if [ -f "$SORTED_BAM" ]; then
                                echo "Successfully created sorted BAM: $SORTED_BAM" | tee -a "$LOG_FILE"
                              else
                                echo "Error: Sorted BAM not created: $SORTED_BAM" | tee -a "$LOG_FILE" >&2
                              fi

                              # Check if the number of unique reads matches N_READS_INPUT
                              if [[ $(samtools view "$SORTED_BAM" | cut -f1 | sort | uniq | wc -l) -ne "$N_READS_INPUT" ]]; then
                                echo "Error: The number of unique reads does not match N_READS_INPUT ($N_READS_INPUT)."
                                exit 1
                              fi
		              done

                             done
                           done
                         done
                       done
                     done
                   done
    #             done
    #           done
    #         done
    #       done
    #     done
    #   done
    # done
  done
done

echo "Script completed successfully."

