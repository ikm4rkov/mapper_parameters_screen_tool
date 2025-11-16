cript to build indexes for different mappers

# Usage: ./build_indexes.sh <genome_file> <gtf_file> <output_prefix> <threads>

# Check if the correct number of arguments is provided
if [ $# -ne 4 ]; then
  echo "Usage: $0 <genome_file> <gtf_file> <output_prefix> <threads>"
  exit 1
fi

genome_file="$1"
gtf_file="$2"
output_prefix="$3"
threads="$4"

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


echo "Starting index building..."

# 1. Bowtie2 Index
echo "Building Bowtie2 index..."
bowtie2-build --threads "$threads" "$genome_file" "$output_prefix.bt2"
if [ $? -ne 0 ]; then
  echo "Error: Bowtie2 index build failed."
  exit 1
fi

# 2. STAR Index
echo "Building STAR index..."
STAR --runThreadN "$threads" --runMode genomeGenerate --genomeDir "$output_prefix.star" --genomeFastaFiles "$genome_file" --sjdbGTFfile "$gtf_file"
if [ $? -ne 0 ]; then
  echo "Error: STAR index build failed."
  exit 1
fi

# 3. HISAT2 Index
echo "Building HISAT2 index..."

# Create exon map file
hisat2_extract_exons.py "$gtf_file" > exons.txt
if [ $? -ne 0 ]; then
    echo "Error: hisat2_extract_exons.py failed."
    exit 1
fi

# Create splicing site map file
hisat2_extract_splice_sites.py "$gtf_file" > splice_sites.txt
if [ $? -ne 0 ]; then
    echo "Error: hisat2_extract_splice_sites.py failed."
    exit 1
fi

# Build HISAT2 index
hisat2-build -p "$threads" --ss splice_sites.txt --exon exons.txt "$genome_file" "$output_prefix.h2"
rm exons.txt splice_sites.txt # Clean up temporary files
if [ $? -ne 0 ]; then
  echo "Error: HISAT2 index build failed."
  exit 1
fi

# 4. BWA-MEM2 Index
echo "Building BWA-MEM2 index..."
bwa-mem2 index -p "$output_prefix.bwa" "$genome_file"
if [ $? -ne 0 ]; then
  echo "Error: BWA-MEM2 index build failed."
  exit 1
fi

echo "Index building complete!"

