#Clean and assemble transcriptomes using raw reads from Illumina RNA bulk sequencing

#Step 1: Create an environment and install the required packages

conda create -n cut
conda activate cut
conda install bioconda::cutadapt
conda install bioconda::trim-galore

#Step 2: Clean paired-end raw reads using TrimGalore

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate cut

/path/to/trim_galore --paired \
--retain_unpaired \
--phred33 \
--output_dir /path/to/output/directory \
--length 36 \
-q 5 \
--stringency 1 \
-e 0.1 \
--fastqc /path/to/forward_raw_reads.fastq /path/to/reverse_raw_reads.fastq #forward and reverse raw reads

# Deactivate Conda environment
conda deactivate

#Step 3: Perform transcriptome assembly using Trinity

conda create -n trin
conda activate trin
conda install bioconda::trinity

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

/path/to/Trinity \
  --seqType fq \
  --left /path/to/forward_clean_reads.fq \
  --right /path/to/reverse_clean_reads.fq \
  --CPU 66 \
  --max_memory 20G \
  --output /path/to/output/directory/your_species_trinity_out_dir

# Deactivate Conda environment
conda deactivate


#Step 4 (Optional): If multiple replicates of forward and reverse reads need to be assembled with Trinity, the previous method may not work. An alternative script is provided

#Step 4.1: Merged replicates for the same genotype

cat forward_replicate_1.fastq forward_replicate_2.fastq > merged_forward_reads.fastq

cat reverse_replicate_1.fastq reverse_replicate_2.fastq > merged_reverse_reads.fastq

#Step 4.2: Run TrimGalore with merged reads as described above

#Step 4.3: Normalize reads

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=900:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

/path/to/insilico_read_normalization.pl \
--seqType fq \
--JM 100G \
--max_cov 30 \
--left /path/to/merged_forward_clean_reads.fq \
--right /path/to/merged_reverse_clean_reads.fq \
--pairs_together \
--PARALLEL_STATS \
--CPU 10 \
--output /path/to/output/directory/normalized_reads

# Deactivate Conda environment
conda deactivate

#Step 4.4: Run Trinity assembly

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

/path/to/Trinity \
  --seqType fq \
  --left /path/to/normalized_reads/your_species_forward_normalized_reads.fq \
  --right /path/to/normalized_reads/your_species_reverse_normalized_reads.fq \
  --CPU 66 \
  --max_memory 20G \
  --output /path/to/output/directory/your_species_trinity_out_dir

# Deactivate Conda environment
conda deactivate