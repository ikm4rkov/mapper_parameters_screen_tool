#script to build indexes for different mappers


# Usage: ./indexing.sh <genome_file> <gtf_file> <threads> <config_file>
# Config file is required as the 4th argument to the script; it must define WORK_DIR and OUTPUT_PREFIX which will be used
# to create WORK_DIR/genomes/{bowtie2,star,hisat2,bwa} where indexes are written.
# Example: ./indexing.sh "/mnt/scratch/rnachrom/ryabykh2018/grid_pig/genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa" "/mnt/scratch/rnachrom/ryabykh2018/grid_pig/genes_GTF/Sus_scrofa.Sscrofa11.1.113.chr.gtf" 8 /mnt/scratch/rnachrom/ryabykh2018/grid_pig/sus_scrofa_mapper_parameters_screening/mappers.conf

# Check if the correct number of arguments is provided (config_file is required)
if [ $# -ne 4 ]; then
  echo "Usage: $0 <genome_file> <gtf_file> <threads> <config_file>"
  exit 1
fi

genome_file="$1"
gtf_file="$2"
threads="$3"
config_file="$4"

# Check if the genome file exists
if [ ! -f "$genome_file" ]; then
  echo "Error: Genome file '$genome_file' not found."
  exit 1
fi

# Check if the GTF file exists
if [ ! -f "$gtf_file" ]; then
  echo "Error: GTF file '$gtf_file' not found."
  exit 1
fi

# Check config file existence and source it (required)
if [ ! -f "$config_file" ]; then
  echo "Error: Config file '$config_file' not found. This parameter is required."
  exit 1
fi
# shellcheck source=/dev/null
source "$config_file"


echo "Starting index building..."


# WORK_DIR and OUTPUT_PREFIX must be defined in the config file
if [ -z "${WORK_DIR:-}" ]; then
  echo "Error: WORK_DIR not set in config file. Please set export WORK_DIR=\"/path/to/work_dir\" in the config."
  exit 1
fi
if [ -z "${OUTPUT_PREFIX:-}" ]; then
  echo "Error: OUTPUT_PREFIX not set in config file. Please set export OUTPUT_PREFIX=\"prefix\" in the config."
  exit 1
fi
work_dir="$WORK_DIR"
output_prefix="$OUTPUT_PREFIX"

# Create work directories for indexes
GENOMES_DIR="$work_dir/genomes"
mkdir -p "$GENOMES_DIR/bowtie2" "$GENOMES_DIR/star" "$GENOMES_DIR/hisat2" "$GENOMES_DIR/bwa"


# 1. Bowtie2 Index
echo "Building Bowtie2 index..."
# Expect BOWTIE2 in config to be a directory containing bowtie2 binaries; construct build command by prefixing
if [ -z "$BOWTIE2" ]; then
  echo "Error: BOWTIE2 not set in config file. Expected a directory path."
  exit 1
fi

BOWTIE2_BUILD_CMD="$BOWTIE2/bowtie2-build"

"$BOWTIE2_BUILD_CMD" --threads "$threads" "$genome_file" "$GENOMES_DIR/bowtie2/$output_prefix"
if [ $? -ne 0 ]; then
  echo "Error: Bowtie2 index build failed."
  exit 1
fi

# 2. STAR Index
echo "Building STAR index..."
# Expect STAR in config to be a directory containing the STAR binary
if [ -z "$STAR" ]; then
  echo "Error: STAR not set in config file. Expected a directory path."
  exit 1
fi
STAR_CMD="$STAR/STAR"

"$STAR_CMD" --runThreadN "$threads" --runMode genomeGenerate --genomeDir "$GENOMES_DIR/star" --genomeFastaFiles "$genome_file" --sjdbGTFfile "$gtf_file"
if [ $? -ne 0 ]; then
  echo "Error: STAR index build failed."
  exit 1
fi

# 3. HISAT2 Index
echo "Building HISAT2 index..."
# Expect HISAT2 in config to be a directory containing hisat2 binaries and helper scripts
if [ -z "$HISAT2" ]; then
  echo "Error: HISAT2 not set in config file. Expected a directory path."
  exit 1
fi

HISAT2_DIR="$HISAT2"
HISAT2_EXTRACT_EXONS="${HISAT2_EXTRACT_EXONS:-$HISAT2_DIR/hisat2_extract_exons.py}"
HISAT2_EXTRACT_SPLICE="${HISAT2_EXTRACT_SPLICE:-$HISAT2_DIR/hisat2_extract_splice_sites.py}"
HISAT2_BUILD_CMD="${HISAT2_BUILD:-$HISAT2_DIR/hisat2-build}"

# create exon map
"$HISAT2_EXTRACT_EXONS" "$gtf_file" > $GENOMES_DIR/hisat2/exons.txt
if [ $? -ne 0 ]; then
    echo "Error: hisat2_extract_exons.py failed."
    exit 1
fi

# create splice sites map
"$HISAT2_EXTRACT_SPLICE" "$gtf_file" > $GENOMES_DIR/hisat2/splice_sites.txt
if [ $? -ne 0 ]; then
    echo "Error: hisat2_extract_splice_sites.py failed."
    exit 1
fi

# Build HISAT2 index (place files into work_dir)
"$HISAT2_BUILD_CMD" -p "$threads" --ss $GENOMES_DIR/hisat2/splice_sites.txt --exon $GENOMES_DIR/hisat2/exons.txt "$genome_file" "$GENOMES_DIR/hisat2/$output_prefix"
rm exons.txt splice_sites.txt # Clean up temporary files
if [ $? -ne 0 ]; then
  echo "Error: HISAT2 index build failed."
  exit 1
fi

# 4. BWA-MEM Index
echo "Building BWA-MEM index..."
# Expect BWA_MEM in config to be a directory containing bwa-mem
if [ -z "$BWA_MEM" ]; then
  echo "Error: BWA_MEM not set in config file. Expected a directory path."
  exit 1
fi
BWA_MEM2_CMD="$BWA_MEM/bwa"
"$BWA_MEM2_CMD" index -p "$GENOMES_DIR/bwa/$output_prefix" "$genome_file"
if [ $? -ne 0 ]; then
  echo "Error: BWA-MEM2 index build failed."
  exit 1
fi

echo "Index building complete!"

