#!/bin/bash

# =========================
# Configuration
# =========================
CONTACT_SCRIPT="/gpfs/ryabykhgrigory/bin/Bam_to_Contacts_prerelease2.py"
MODE="STAR"
ENZYME="ATA"
PAIRS_FILE="pairs.tsv"
RUN_MODE="all"   # default mode

# =========================
# Argument parsing
# =========================
print_usage() {
    echo "Usage: $0 [--mode all|one] [--pairs FILE]"
    echo
    echo "  --mode   all | one    : 'all' = full cross of all RNA/DNA BAMs (default)"
    echo "                          'one' = only RNA/DNA pairs in each line of pairs.tsv"
    echo "  --pairs  FILE         : path to a custom pairs.tsv file"
    echo
    exit 1
}

# Parse arguments with getopt
TEMP=$(getopt -o "" --long mode:,pairs:,help -n "$0" -- "$@")
if [[ $? != 0 ]]; then print_usage; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        --mode)
            RUN_MODE="$2"
            shift 2
            ;;
        --pairs)
            PAIRS_FILE="$2"
            shift 2
            ;;
        --help)
            print_usage
            ;;
        --)
            shift
            break
            ;;
        *)
            print_usage
            ;;
    esac
done

# =========================
# Validate input
# =========================
if [[ ! -f "$PAIRS_FILE" ]]; then
    echo "❌ Error: pairs file not found: $PAIRS_FILE"
    exit 1
fi

if [[ "$RUN_MODE" != "all" && "$RUN_MODE" != "one" ]]; then
    echo "❌ Error: invalid mode '$RUN_MODE' (expected 'all' or 'one')."
    exit 1
fi

echo "🧩 Mode: $RUN_MODE"
echo "📄 Pairs file: $PAIRS_FILE"
echo

# =========================
# Extract BAM paths
# =========================
rna_files=($(cut -f1 "$PAIRS_FILE" | sed 's/_log.txt$/_alignments_sorted.bam/' | sort -u))
dna_files=($(cut -f2 "$PAIRS_FILE" | sed 's/_log.txt$/_alignments_sorted.bam/' | sort -u))

echo "📄 Unique RNA BAM files: ${#rna_files[@]}"
echo "📄 Unique DNA BAM files: ${#dna_files[@]}"
echo

# =========================
# Process pairs
# =========================

if [[ "$RUN_MODE" == "all" ]]; then
    # Full cross product of RNA and DNA files
    for rna_file in "${rna_files[@]}"; do
        if [[ ! -f "$rna_file" ]]; then
            echo "⚠️  RNA BAM not found: $rna_file"
            continue
        fi
        for dna_file in "${dna_files[@]}"; do
            if [[ ! -f "$dna_file" ]]; then
                echo "⚠️  DNA BAM not found: $dna_file"
                continue
            fi

            rna_base=$(basename "$rna_file" _alignments_sorted.bam)
            dna_base=$(basename "$dna_file" _alignments_sorted.bam)
            rna_dir=$(basename "$(dirname "$rna_file")")
            dna_dir=$(basename "$(dirname "$dna_file")")
            output_prefix="rna_${rna_dir}_${rna_base}__dna_${dna_dir}_${dna_base}"

            echo "🔄 Intersecting:"
            echo "  RNA: $rna_file"
            echo "  DNA: $dna_file"
            echo "  Output prefix: $output_prefix"

            python "$CONTACT_SCRIPT" \
                -r1 "$rna_file" \
                -r2 "$dna_file" \
                -mr "$MODE" \
                -md "$MODE" \
                -e "$ENZYME" \
                -p "$output_prefix"
        done
    done

else
    # One-to-one mode: process pairs line-by-line
    while IFS=$'\t' read -r rna_log dna_log; do
        rna_file="${rna_log/_log.txt/_alignments_sorted.bam}"
        dna_file="${dna_log/_log.txt/_alignments_sorted.bam}"

        if [[ ! -f "$rna_file" ]]; then
            echo "⚠️  RNA BAM not found: $rna_file"
            continue
        fi
        if [[ ! -f "$dna_file" ]]; then
            echo "⚠️  DNA BAM not found: $dna_file"
            continue
        fi

        rna_base=$(basename "$rna_file" _alignments_sorted.bam)
        dna_base=$(basename "$dna_file" _alignments_sorted.bam)
        rna_dir=$(basename "$(dirname "$rna_file")")
        dna_dir=$(basename "$(dirname "$dna_file")")
        output_prefix="rna_${rna_dir}_${rna_base}__dna_${dna_dir}_${dna_base}"

        echo "🔄 Intersecting:"
        echo "  RNA: $rna_file"
        echo "  DNA: $dna_file"
        echo "  Output prefix: $output_prefix"

        python "$CONTACT_SCRIPT" \
            -r1 "$rna_file" \
            -r2 "$dna_file" \
            -mr "$MODE" \
            -md "$MODE" \
            -e "$ENZYME" \
            -p "$output_prefix"

    done < "$PAIRS_FILE"
fi

