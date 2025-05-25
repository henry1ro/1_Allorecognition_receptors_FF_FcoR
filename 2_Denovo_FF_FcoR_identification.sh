# De novo isolation of Fester and Fester coreceptor genes from a transcriptome

# Since there are over 100 allorecognition genes (fuhc, Fester, and Fester coreceptor),
# manually extracting sequences that BLAST-match each gene is impractical.
# The following pipeline automates the process.

# Step 1: Perform tBLASTx with allorecognition genes

makeblastdb -in your_transcriptome.fasta -dbtype nucl -out your_transcriptome_name

tblastx -query allorecognition_genes.fasta -db your_transcriptome_name -out blast_output.xml -outfmt 5

# Step 2: Extract headers for transcripts with E-value below the threshold

# Create and activate a conda environment with required tools
conda create -n blast
conda activate blast
conda install -c bioconda biopython seqtk seqkit transdecoder

# Save this file as blast_script.py
from Bio.Blast import NCBIXML

with open("blast_output.xml", "r") as f:  # Replace with your BLAST output
    blast_records = NCBIXML.parse(f)
    with open("out_your_transcriptome.txt", "w") as out_file:  # Output file for transcript headers
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.1:  # Adjust threshold as needed (0.1 or 0.01)
                        out_file.write(alignment.title.split(" ")[0] + "\n")  # Save header only

# Run this script:
# python blast_script.py

# Step 3: Extract IDs from the text file
# Use a spreadsheet application (e.g., Excel) or command-line tools to extract IDs and save as a text file.

# Step 4: Extract sequences
seqtk subseq /path/to/your_transcriptome.fasta /path/to/your_IDs_to_extract.txt > /path/to/extracted_sequences.fasta

# Step 5: Remove duplicate sequences
seqkit rmdup -s /path/to/extracted_sequences.fasta > /path/to/extracted_sequences_non_redundant.fasta

# Step 6: Predict ORFs
TransDecoder.LongOrfs -t /path/to/extracted_sequences_non_redundant.fasta

# Step 7 (Optional): Remove duplicate sequences
seqkit rmdup -s /path/to/longest_orfs.cds > /path/to/longest_orfs_non_redundant.cds

# Step 8: Predict protein-coding regions
TransDecoder.Predict --cpu 20 --single_best_only -t /path/to/extracted_sequences_non_redundant.fasta

# Step 9: Remove duplicate sequences
seqkit rmdup -s /path/to/your_sequences.transdecoder.cds > /path/to/your_sequences.transdecoder_non_redundant.cds
