#!/bin/bash
# Usage: ./script.sh dummy_arg /path/to/star_montecarlo_outputs

PATH_TO_BASE="$1"

if [ -z "$PATH_TO_BASE" ]; then
  echo "Error: missing path argument."
  echo "Usage: $0 dummy_arg /path/to/base_directory"
  exit 1
fi

# Try to find existing BAM files first
bam_files=$(find "$PATH_TO_BASE" -type f -name "*.bam")

if [ -n "$bam_files" ]; then
  echo "Found existing BAM files — sorting by read ID..."
  echo "$bam_files" | parallel -j 5 '
    samtools sort -n {} -o {.}_sorted.bam
  '
else
  echo "No BAM files found — converting SAM → BAM, then sorting..."
  find "$PATH_TO_BASE" -type f -name "*.sam" | parallel -j 5 '
    samtools view -b {} > {.}_tmp.bam &&
    samtools sort -n {.}_tmp.bam -o {.}_sorted.bam &&
    rm {.}_tmp.bam
  '
fi

echo "Processing complete."

