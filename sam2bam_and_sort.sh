#!/bin/bash
# Usage: ./script.sh dummy_arg /path/to/star_montecarlo_outputs

PATH_TO_BASE="$1"

if [ -z "$PATH_TO_BASE" ]; then
  echo "Error: missing path argument."
  echo "Usage: $0 dummy_arg /path/to/base_directory"
  exit 1
fi

# Find BAM files that DO NOT contain "_sorted" in their filename
bam_files=$(find "$PATH_TO_BASE" -type f -name "*.bam" ! -name "*_sorted*.bam")

if [ -n "$bam_files" ]; then
  echo "Found existing BAM files — sorting by read ID (excluding already sorted)..."
  echo "$bam_files" | parallel -j 5 '
    samtools sort -n {} -o {.}_sorted.bam
  '
else
  echo "No unsorted BAM files found — converting SAM → BAM, then sorting (excluding already sorted)..."
  find "$PATH_TO_BASE" -type f -name "*.sam" ! -name "*_sorted*.sam" | parallel -j 5 '
    samtools view -b {} > {.}_tmp.bam &&
    samtools sort -n {.}_tmp.bam -o {.}_sorted.bam &&
    rm {.}_tmp.bam
  '
fi

echo "Processing complete."

