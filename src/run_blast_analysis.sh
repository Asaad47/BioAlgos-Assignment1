#!/bin/bash

# Set paths for reference genomes
GENOME_FILES=(
    "../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
    "../data/2_bsub_ncbi_dataset/ncbi_dataset/data/GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna"
    "../data/3_paer_ncbi_dataset/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
    "../data/4_saur_ncbi_dataset/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna"
    "../data/5_mtub_ncbi_dataset/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna"
)

# Concatenate all genome files into one reference file
echo "Concatenating genome files into five_genomes.fna..."
cat "${GENOME_FILES[@]}" > ../data/five_genomes.fna

# Build BLAST database
echo "Building BLAST database..."
makeblastdb -in ../data/five_genomes.fna -dbtype nucl -out reference_db

# Set paths for sequencing reads
READ_FILES=(
    "../data/sequence_reads/simulated_reads_no_errors_10k_R1.fastq"
    "../data/sequence_reads/simulated_reads_no_errors_10k_R2.fastq"
    "../data/sequence_reads/simulated_reads_miseq_10k_R1.fastq"
    "../data/sequence_reads/simulated_reads_miseq_10k_R2.fastq"
)

# Convert FASTQ files to FASTA (if not already in FASTA format)
echo "Converting FASTQ files to FASTA format..."
for READ_FILE in "${READ_FILES[@]}"; do
    FASTA_FILE="${READ_FILE%.fastq}.fasta"  # Change extension to .fasta
    seqtk seq -A "$READ_FILE" > "$FASTA_FILE"
done

# Run BLAST on each read file
echo "Running BLAST for each read file..."
for READ_FILE in "${READ_FILES[@]}"; do
    FASTA_FILE="${READ_FILE%.fastq}.fasta"
    OUTPUT_FILE="${FASTA_FILE%.fasta}_blast_results.txt"
    if [ ! -f "$OUTPUT_FILE" ]; then
        echo "Output file $OUTPUT_FILE not found. Creating empty file..."
        touch "$OUTPUT_FILE"
    fi

    echo "Processing $FASTA_FILE..."
    /usr/bin/time -l blastn -query "$FASTA_FILE" -db reference_db -out "$OUTPUT_FILE" -outfmt 6 -max_target_seqs 5 -evalue 1e-5
done

# Generate summary report
echo "Generating classification summary..."
for READ_FILE in "${READ_FILES[@]}"; do
    FASTA_FILE="${READ_FILE%.fastq}.fasta"
    OUTPUT_FILE="${FASTA_FILE%.fasta}_blast_results.txt"

    echo "Summary for $OUTPUT_FILE:"
    cut -f2 "$OUTPUT_FILE" | sort | uniq -c
    echo "------------------------------------------------"
done

echo "Analysis complete!"
