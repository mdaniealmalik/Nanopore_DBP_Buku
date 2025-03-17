#!/bin/bash

# Directory where FASTQ files are located
input_dir="1_Sample"

# Output directory where filtered files will be saved
output_dir="2_NanoFilt_output"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Quality threshold for NanoFilt
quality_threshold=10

# Minimum length for NanoFilt
min_length=180

# Maximum length for NanoFilt
max_length=320

# Loop through each FASTQ file in the input directory
for fastq_file in "$input_dir"/*.fastq.gz; do
    # Extract the base name of the file (without directory and extension)
    base_name=$(basename "$fastq_file" .fastq.gz)
    
    # Set the output file name
    output_file="$output_dir/${base_name}_filtered.fastq"
    
    # Run NanoFilt to filter the FASTQ file
    gunzip -c "$fastq_file" | NanoFilt -q "$quality_threshold" -l "$min_length" --maxlength "$max_length" | gzip  > "$output_file"
    
    # Print a message indicating the file has been processed
    echo "Processed $fastq_file and saved filtered data to $output_file"
done

# Directory containing the input files
INPUT_DIR="2_NanoFilt_output"  # Change this to your input directory if different

# Directory for output files
OUTPUT_DIR="3_cutadapt_output"  # Change this to your desired output directory

# Cutadapt parameters
PRIMER_SEQUENCE_FORWARD="TTTCTGTTGGTGCTGATATTGCGTTGGTAAATCTCGTGCCAGC"
PRIMER_SEQUENCE_REVERSE="ACTTGCCTGTCGCTCTATCTTCCATAGTGGGGTATCTAATCCTAGTTTG"
ERROR_RATE="0.2"        # Set this to the desired error rate
MIN_LENGTH="150"         # Set this to the minimum length of reads to keep
MAX_LENGTH="200"        # Set this to the maximum length of reads to keep

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over all FASTQ files in the input directory
for FILE in "$INPUT_DIR"/*.fastq; do
  # Get the base name of the file without the directory
  BASENAME=$(basename "$FILE" .fastq)
  
  # Define the output file name
  OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_trimmed_forward.fastq"

  # Run cutadapt with the specified options
  cutadapt -g "$PRIMER_SEQUENCE_FORWARD" \
           -o "$OUTPUT_FILE" \
           -e "$ERROR_RATE" \
           -m "$MIN_LENGTH" \
           -M "$MAX_LENGTH" \
           --discard-untrimmed \
           "$FILE"

  # Print status message
  echo "Processed $FILE -> $OUTPUT_FILE"
done


# Loop over all FASTQ files in the input directory
for FILE in "$INPUT_DIR"/*.fastq; do
  # Get the base name of the file without the directory
  BASENAME=$(basename "$FILE" .fastq)
  
  # Define the output file name
  OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_trimmed_reverse.fastq"

  # Run cutadapt with the specified options
  cutadapt -g "$PRIMER_SEQUENCE_REVERSE" \
           -o "$OUTPUT_FILE" \
           -e "$ERROR_RATE" \
           -m "$MIN_LENGTH" \
           -M "$MAX_LENGTH" \
           --discard-untrimmed \
           "$FILE"

  # Print status message
  echo "Processed $FILE -> $OUTPUT_FILE"
done


# Define the input and output directories
input_dir="3_cutadapt_output"
output_dir="4_combined_fastq_from_cutadapt"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over all forward FASTQ files in the input directory
for forward_file in ${input_dir}/*_forward.fastq; do
  # Extract the base name (excluding "_forward.fastq")
  base_name=$(basename "$forward_file" "_forward.fastq")

  # Define the corresponding reverse file and the output combined file
  reverse_file="${input_dir}/${base_name}_reverse.fastq"
  combined_file="${output_dir}/${base_name}_combined.fastq"

  # Check if the corresponding reverse file exists
  if [[ -f "$reverse_file" ]]; then
    echo "Combining $forward_file and $reverse_file into $combined_file"
    
    # Combine the files (concatenation in this case)
    cat "$forward_file" "$reverse_file" > "$combined_file"
  else
    echo "Reverse file $reverse_file not found for $forward_file"
  fi
done

echo "Combination process completed. Combined files are in $output_dir."

# Define the input and output directories
input_dir="4_combined_fastq_from_cutadapt"
output_dir="5_vsearch"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over all FASTQ files in the input directory
for fastq_file in ${input_dir}/*.fastq; do
  # Extract the base name (excluding the file extension)
  base_name=$(basename "$fastq_file" ".fastq")

  # Define the output FASTA file name
  fasta_file="${output_dir}/${base_name}.fasta"

  # Convert FASTQ to FASTA using seqtk
  echo "Converting $fastq_file to $fasta_file"
  seqtk seq -A "$fastq_file" > "$fasta_file"
done

echo "Conversion process completed. FASTA files are in $output_dir."

# Directory containing FASTA files
input_dir="5_vsearch"

# Output directory for renamed FASTA files
output_dir="5_vsearch/rename_fasta"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over each FASTA file in the directory
for fasta_file in "$input_dir"/*.fasta; do
  # Check if the file exists
  if [[ -f "$fasta_file" ]]; then
    # Extract the base name of the file (without path and extension)
    base_name=$(basename "$fasta_file" .fasta)
    
    # Output file path
    output_file="$output_dir/${base_name}_rename.fasta"

    # Count the number of sequences
    sequence_count=$(sed '/^>/d' "$fasta_file" | wc -l)

    # Generate new headers and rename them
    {
      for i in $(seq 1 "$sequence_count"); do
        echo "$base_name;$i"
      done
    } | paste - <(sed '/^>/d' "$fasta_file") | sed -e 's/^/>/' -e 's/\t/\n/' > "$output_file"

    echo "Headers for $fasta_file have been renamed and saved to $output_file"
  else
    echo "File $fasta_file does not exist."
  fi
done

#file for combine all fasta
mkdir -p 5_vsearch/rename_fasta/combine

#combine all barcode fasta into single fasta
cat 5_vsearch/rename_fasta/* >> 5_vsearch/rename_fasta/combine/combine.fasta

#dereplicate
vsearch --derep_fulllength 5_vsearch/rename_fasta/combine/combine.fasta --output 5_vsearch/rename_fasta/combine/dereplicated_combine.fasta --sizeout

#cluster
vsearch --cluster_fast 5_vsearch/rename_fasta/combine/dereplicated_combine.fasta --centroids 5_vsearch/rename_fasta/combine/cluster_combine.fasta --id 0.95 --sizein --sizeout

#chimera detection
vsearch --uchime_denovo 5_vsearch/rename_fasta/combine/cluster_combine.fasta --chimeras 5_vsearch/rename_fasta/combine/chimeras_combine.fasta --nonchimeras 5_vsearch/rename_fasta/combine/nonchimeras_combine.fasta

# Process the file: remove semicolons from headers
awk '/^>/ { gsub(/;/, "", $0) } { print }' 5_vsearch/rename_fasta/combine/nonchimeras_combine.fasta > 5_vsearch/rename_fasta/combine/nonchimeras_combine_rename.fasta

# otu table
vsearch --usearch_global 5_vsearch/rename_fasta/combine/combine.fasta  \
    --db 5_vsearch/rename_fasta/combine/nonchimeras_combine_rename.fasta --id 0.95 --threads 4 \
    --otutabout otu_table.tsv --sizein --strand plus

# blastn
makeblastdb -in database/database.fasta -dbtype nucl -out database/database

blastn -query 5_vsearch/rename_fasta/combine/nonchimeras_combine_rename.fasta -db database/database -evalue '0.001' -out result_blastn.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 10 -strand both -dust yes -max_target_seqs 2 -perc_identity 90 -qcov_hsp_perc 90


# building taxon_table
# Get the current working directory
WD=$(pwd)
# Run R script with the working directory argument
Rscript -e "setwd('$WD'); source('script_r_taxon_table.R')"
