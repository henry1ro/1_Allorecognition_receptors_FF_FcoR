#Prepare the reference sequences

dos2unix /path/to/3_LONG_all_genes_splicing_version_2.fasta

#Mapping the raw reads to full-lenght and alternative splicing variants of FF genes

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

align_and_estimate_abundance.pl \
--est_method RSEM \
--aln_method bowtie2 \
--samples_file /path/to/SampleFile_TotalBotryllus_6.txt \
--transcripts /path/to/3_LONG_all_genes_splicing_version_2.fasta \
--seqType fq \
--prep_reference \
--output_dir /path/to/output/directory

# Deactivate Conda environment
conda deactivate

#Prepare the matrix

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

abundance_estimates_to_matrix.pl \
--est_method RSEM \
--out_prefix all \
--gene_trans_map none \
--name_sample_by_basedir \
/path/2y_SB825_2/RSEM.genes.results \
/path/2y_SB825_3/RSEM.genes.results \
/path/2y_SMGO_1/RSEM.genes.results \
/path/2y_SMGO_2/RSEM.genes.results \
/path/2y_SMGO_3/RSEM.genes.results \
/path/2y_SMH_1/RSEM.genes.results \
/path/2y_SMH_2/RSEM.genes.results \
/path/2y_SMH_3/RSEM.genes.results \
/path/3y_SB825_1/RSEM.genes.results \
/path/3y_SB825_2/RSEM.genes.results \
/path/3y_SB825_3/RSEM.genes.results \
/path/3y_SMGO_1/RSEM.genes.results \
/path/3y_SMGO_2/RSEM.genes.results \
/path/3y_SMGO_3/RSEM.genes.results \
/path/3y_SMH_1/RSEM.genes.results \
/path/3y_SMH_2/RSEM.genes.results \
/path/3y_SMH_3/RSEM.genes.results \
/path/4m_SB825_1/RSEM.genes.results \
/path/4m_SB825_2/RSEM.genes.results \
/path/4m_SB825_3/RSEM.genes.results \
/path/4m_SMGO_1/RSEM.genes.results \
/path/4m_SMGO_2/RSEM.genes.results \
/path/4m_SMGO_3/RSEM.genes.results \
/path/4m_SMH_1/RSEM.genes.results \
/path/4m_SMH_2/RSEM.genes.results \
/path/4m_SMH_3/RSEM.genes.results \
/path/A1_Fer_801/RSEM.genes.results \
/path/A1_Fer_841/RSEM.genes.results \
/path/A1_Inf_001/RSEM.genes.results \
/path/A1_Inf_802/RSEM.genes.results \
/path/A1_MixFer/RSEM.genes.results \
/path/A1_MixInf/RSEM.genes.results \
/path/A2_Fer_801/RSEM.genes.results \
/path/A2_Fer_841/RSEM.genes.results \
/path/A2_Inf_001/RSEM.genes.results \
/path/A2_Inf_802/RSEM.genes.results \
/path/A2_MixFer/RSEM.genes.results \
/path/A2_MixInf/RSEM.genes.results \
/path/B1_Fer_801/RSEM.genes.results \
/path/B1_Fer_841/RSEM.genes.results \
/path/B1_Inf_001/RSEM.genes.results \
/path/B1_Inf_802/RSEM.genes.results \
/path/B1_MixFer/RSEM.genes.results \
/path/B1_MixInf/RSEM.genes.results \
/path/B2_Fer_801/RSEM.genes.results \
/path/B2_Fer_841/RSEM.genes.results \
/path/B2_Inf_001/RSEM.genes.results \
/path/B2_Inf_802/RSEM.genes.results \
/path/B2_MixFer/RSEM.genes.results \
/path/B2_MixInf/RSEM.genes.results \
/path/C1_Fer_801/RSEM.genes.results \
/path/C1_Fer_841/RSEM.genes.results \
/path/C1_Inf_001/RSEM.genes.results \
/path/C1_Inf_802/RSEM.genes.results \
/path/C1_MixFer/RSEM.genes.results \
/path/C1_MixInf/RSEM.genes.results \
/path/C2_Fer_801/RSEM.genes.results \
/path/C2_Fer_841/RSEM.genes.results \
/path/C2_Inf_001/RSEM.genes.results \
/path/C2_Inf_802/RSEM.genes.results \
/path/C2_MixFer/RSEM.genes.results \
/path/C2_MixInf/RSEM.genes.results \
/path/D_Fer_801/RSEM.genes.results \
/path/D_Fer_841/RSEM.genes.results \
/path/D_Inf_001/RSEM.genes.results \
/path/D_Inf_802/RSEM.genes.results \
/path/D_MixFer/RSEM.genes.results \
/path/D_MixInf/RSEM.genes.results \
/path/ALDH_Loser_1_minus/RSEM.genes.results \
/path/ALDH_Loser_1_plus/RSEM.genes.results \
/path/ALDH_Loser_7_minus/RSEM.genes.results \
/path/ALDH_Loser_7_plus/RSEM.genes.results \
/path/ALDH_Winner_2B_minus/RSEM.genes.results \
/path/ALDH_Winner_2B_plus/RSEM.genes.results \
/path/ALDH_Winner_6_minus/RSEM.genes.results \
/path/ALDH_Winner_6_plus/RSEM.genes.results \
/path/Ampullae_SB825FDA_Rep1/RSEM.genes.results \
/path/Ampullae_SB825FDA_Rep2/RSEM.genes.results \
/path/Ampullae_SB825FDA_Rep3/RSEM.genes.results \
/path/Ampullae_SB825FEA_Rep1/RSEM.genes.results \
/path/Ampullae_SB825FEA_Rep2/RSEM.genes.results \
/path/Ampullae_SB825FEA_Rep3/RSEM.genes.results \
/path/Ampullae_SB825FX_Rep1/RSEM.genes.results \
/path/Ampullae_SB825FX_Rep2/RSEM.genes.results \
/path/Ampullae_SB825FX_Rep3/RSEM.genes.results \
/path/Ampullae_SMGO_Rep2/RSEM.genes.results \
/path/Ampullae_SMGO_Rep3/RSEM.genes.results \
/path/Ampullae_SMGO_Technical_Rep1/RSEM.genes.results \
/path/Ampullae_SMGO_Technical_Rep2/RSEM.genes.results \
/path/Ampullae_SMH_Rep1/RSEM.genes.results \
/path/Ampullae_SMH_Rep2/RSEM.genes.results \
/path/Ampullae_SMH_Rep3/RSEM.genes.results \
/path/Blood_FL3/RSEM.genes.results \
/path/Blood_FLLC2/RSEM.genes.results \
/path/Blood_LFC4/RSEM.genes.results \
/path/IA6_79C_plus_rep1/RSEM.genes.results \
/path/IA6_79C_plus_rep2/RSEM.genes.results \
/path/IA6_79G_plus_rep1/RSEM.genes.results \
/path/IA6_79G_plus_rep2/RSEM.genes.results \
/path/IA6_81A_minu_rep1/RSEM.genes.results \
/path/IA6_81A_minu_rep2/RSEM.genes.results \
/path/IA6_81A_plus_rep1/RSEM.genes.results \
/path/IA6_81A_plus_rep2/RSEM.genes.results \
/path/IA6_81C_minu_rep1/RSEM.genes.results \
/path/IA6_81C_minu_rep2/RSEM.genes.results \
/path/IA6_81C_plus_rep1/RSEM.genes.results \
/path/IA6_81C_plus_rep2/RSEM.genes.results \
/path/Male_M001_A1/RSEM.genes.results \
/path/Male_M001_A2/RSEM.genes.results \
/path/Male_M001_B1/RSEM.genes.results \
/path/Male_M001_B2/RSEM.genes.results \
/path/Male_M001_C1/RSEM.genes.results \
/path/Male_M001_C2/RSEM.genes.results \
/path/Male_M001_D/RSEM.genes.results \
/path/Male_MADM_A1/RSEM.genes.results \
/path/Male_MADM_A2/RSEM.genes.results \
/path/Male_MADM_B1/RSEM.genes.results \
/path/Male_MADM_B2/RSEM.genes.results \
/path/Male_MADM_C1/RSEM.genes.results \
/path/Male_MADM_C2/RSEM.genes.results \
/path/Male_MADM_D/RSEM.genes.results \
/path/Male_MDR_A1/RSEM.genes.results \
/path/Male_MDR_A2/RSEM.genes.results \
/path/Male_MDR_B1/RSEM.genes.results \
/path/Male_MDR_B2/RSEM.genes.results \
/path/Male_MDR_C1/RSEM.genes.results \
/path/Male_MDR_C2/RSEM.genes.results \
/path/Male_MDR_D/RSEM.genes.results \
/path/Oocyte_with_follicles/RSEM.genes.results \
/path/Oocyte_without_follicles/RSEM.genes.results \
/path/WB_825_rep1/RSEM.genes.results \
/path/WB_825_rep2/RSEM.genes.results \
/path/WB_825_rep3/RSEM.genes.results \
/path/WB_SKI_rep1/RSEM.genes.results \
/path/WB_SKI_rep2/RSEM.genes.results \
/path/WB_SKI_rep3/RSEM.genes.results \
/path/WB_SMGO_rep1/RSEM.genes.results \
/path/WB_SMGO_rep2/RSEM.genes.results \
/path/WB_SMGO_rep3/RSEM.genes.results \
/path/WB_SMH_802/RSEM.genes.results \
/path/WB_SMH_rep1/RSEM.genes.results \
/path/WB_SMH_rep2/RSEM.genes.results \
/path/T0_rep1/RSEM.genes.results \
/path/T0_rep2/RSEM.genes.results \
/path/T0_rep3/RSEM.genes.results \
/path/T6_rep1/RSEM.genes.results \
/path/T6_rep2/RSEM.genes.results \
/path/T6_rep3/RSEM.genes.results \
/path/T18_rep1/RSEM.genes.results \
/path/T18_rep2/RSEM.genes.results \
/path/T18_rep3/RSEM.genes.results \
/path/T24_rep1/RSEM.genes.results \
/path/T24_rep2/RSEM.genes.results \
/path/T24_rep3/RSEM.genes.results

# Deactivate Conda environment
conda deactivate

#Heatmap

# Step 1: Load the necessary libraries
library(ComplexHeatmap)
library(circlize)  # For color palette customization

# Step 2: Load the data from the text file
data <- read.table("/path/to/all.isoform.TMM.EXPR.txt",    #path to mapping text file results
                   header = TRUE, sep = "\t", row.names = 1)  # Assuming genes are in the first column

# Step 3: Convert the row names (gene names) and ensure data is numeric
gene_names <- rownames(data)  # Extract row names (gene names)
colnames(data) <- gsub("X", "", colnames(data))  # Clean column names if needed (e.g., if they start with X)

# Step 4: Convert the data to numeric (if not already numeric)
data[] <- lapply(data, as.numeric)

# Optional: Replace zero values with a small constant (e.g., 1) for visualization purposes
data[data == 0] <- 1

# Step 5: Optionally, perform log2 transformation of the data
data_log <- log2(data + 1)  # Add 1 to avoid log(0)

# Step 6: Create the heatmap without grouping genes
Heatmap(data_log, 
        name = "Expression",  # Heatmap color legend title
        column_title = "Transcriptomes",  # Title for the columns
        row_title = "Genes",  # Title for the rows
        show_row_names = TRUE,  # Show gene names (row labels)
        show_column_names = TRUE,  # Show transcriptome names (column labels)
        cluster_rows = FALSE,  # Do not cluster rows
        clustering_distance_columns = "euclidean",  # Method for distance calculation in columns
        clustering_method_columns = "complete",  # Method for hierarchical clustering of columns
        column_dend_side = "top",  # Place column dendrogram at the top
        col = colorRamp2(c(0, 5, 10), c("blue", "white", "red")),  # Blue-white-red palette similar to pheatmap
        heatmap_legend_param = list(title = "Expression", at = c(0, 5, 10), labels = c("Low", "Medium", "High")),  # Custom legend
        row_dend_width = unit(5, "cm"),  # Adjust the width of the row dendrogram
        column_dend_height = unit(2, "cm"))  # Reduced the height of the column dendrogram
