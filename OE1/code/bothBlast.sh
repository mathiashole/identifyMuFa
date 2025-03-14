#!/bin/bash

# Verification of arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <query_fasta> <output_dir> <genome_fasta1> [<genome_fasta2> ...]"
    exit 1
fi

# Variables
query_fasta="$1"
output_dir="$2"
genomes="${@:3}"
# output_dir="blast_results"
amino_output="${query_fasta%.fasta}_amino.fasta"

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# 1. Creación de la base de datos BLAST para cada genoma
echo "Creating BLAST databases for genomes..."
for genome in $genomes; do
    genome_name=$(basename "$genome" .fasta)
    makeblastdb -in "$genome" -dbtype nucl -out "$output_dir/$genome_name"
    echo "Database created for $genome"
done

# 2. Búsqueda de nucleótidos contra bases de datos de genomas (blastn)
echo "Running blastn for nucleotide search..."
for genome in $genomes; do
    genome_name=$(basename "$genome" .fasta)
    blastn_out="$output_dir/blastn_${genome_name}.txt"
    
    blastn -query "$query_fasta" -db "$output_dir/$genome_name" -out "$blastn_out" -evalue 1e-5 -outfmt '6 std qlen slen'
    echo "blastn results saved to $blastn_out"
done

# # 2.1 Extraer secuencia de coincidencias
# echo "Running gscissors for nucleotide sequence ..."
# for genome in $genomes; do
#     genome_name=$(basename "$genome" .fasta)
#     blastn_out="$output_dir/blastn_${genome_name}.txt"
#     blastn_out_fasta="$output_dir/blastn_${genome_name}.fasta"
    
#     if [ -s "$blastn_out" ]; then
#         ./gscissors.pl --fasta "$genome" --coordinates "$blastn_out" --format blast --output "$blastn_out_fasta"
#         echo "Fasta nucleotide sequence saved to $blastn_out_fasta"
#     else
#         echo "Warning: No matches found in $blastn_out; skipping extraction."
#     fi
# done

# 3. Conversión de secuencias a aminoácidos (transeq)
echo "Converting nucleotide sequences to amino acids with transeq..."
transeq -clean -sequence "$query_fasta" -outseq "$amino_output"
echo "Amino acid sequences saved to $amino_output"

# 4. Búsqueda de aminoácidos contra bases de datos de genomas (tblastn)
echo "Running tblastn for amino acid search..."
for genome in $genomes; do
    genome_name=$(basename "$genome" .fasta)
    tblastn_out="$output_dir/tblastn_${genome_name}.txt"
    
    tblastn -query "$amino_output" -db "$output_dir/$genome_name" -out "$tblastn_out" -evalue 1e-5 -outfmt '6 std qlen slen'
    echo "tblastn results saved to $tblastn_out"
done

# # 4.1 Extraer secuencia de coincidencias en aminoácidos
# echo "Running gscissors for amino acid sequence ..."
# for genome in $genomes; do
#     genome_name=$(basename "$genome" .fasta)
#     tblastn_out="$output_dir/tblastn_${genome_name}.txt"
#     tblastn_out_fasta="$output_dir/tblastn_${genome_name}.fasta"
    
#     if [ -s "$tblastn_out" ]; then
#         ./gscissors.pl --fasta "$genome" --coordinates "$tblastn_out" --format blast --output "$tblastn_out_fasta"
#         echo "Fasta amino acid sequence saved to $tblastn_out_fasta"
#     else
#         echo "Warning: No matches found in $tblastn_out; skipping extraction."
#     fi
# done

echo "All searches completed. Results saved in the $output_dir directory."