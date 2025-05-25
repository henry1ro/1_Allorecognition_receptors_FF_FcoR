#Mapping raw reads to allorecognition genes with Kallisto

#Step 1: Install Trinity software

conda create -n trin
conda activate trin
conda install bioconda::trinity

#Step 2: Mapping with Kallisto

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

align_and_estimate_abundance.pl \
--est_method kallisto \
--samples_file /path/to/SampleFile_TotalBotryllus_2.txt \    #sample file that specifies the path to the raw reads
--transcripts /path/to/allogenes_References_05_18_FINAL.fasta \  #allorecognition genes reference 
--seqType fq \
--prep_reference \
--output_dir /path/to/output/directory

# Deactivate Conda environment
conda deactivate

#STEP 3: Prepare the matrix

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

abundance_estimates_to_matrix.pl \
--est_method kallisto \
--out_prefix all \
--gene_trans_map none \
--name_sample_by_basedir \
/path/to/2y_SB825_1/abundance.tsv \
/path/to/2y_SB825_2/abundance.tsv \
/path/to/2y_SB825_3/abundance.tsv \
/path/to/2y_SMGO_1/abundance.tsv \
/path/to/2y_SMGO_2/abundance.tsv \
/path/to/2y_SMGO_3/abundance.tsv \
/path/to/2y_SMH_1/abundance.tsv \
/path/to/2y_SMH_2/abundance.tsv \
/path/to/2y_SMH_3/abundance.tsv \
/path/to/3y_SB825_1/abundance.tsv \
/path/to/3y_SB825_2/abundance.tsv \
/path/to/3y_SB825_3/abundance.tsv \
/path/to/3y_SMGO_1/abundance.tsv \
/path/to/3y_SMGO_2/abundance.tsv \
/path/to/3y_SMGO_3/abundance.tsv \
/path/to/3y_SMH_1/abundance.tsv \
/path/to/3y_SMH_2/abundance.tsv \
/path/to/3y_SMH_3/abundance.tsv \
/path/to/4m_SB825_1/abundance.tsv \
/path/to/4m_SB825_2/abundance.tsv \
/path/to/4m_SB825_3/abundance.tsv \
/path/to/4m_SMGO_1/abundance.tsv \
/path/to/4m_SMGO_2/abundance.tsv \
/path/to/4m_SMGO_3/abundance.tsv \
/path/to/4m_SMH_1/abundance.tsv \
/path/to/4m_SMH_2/abundance.tsv \
/path/to/4m_SMH_3/abundance.tsv \
/path/to/A1_Fer_801/abundance.tsv \
/path/to/A1_Fer_841/abundance.tsv \
/path/to/A1_Inf_001/abundance.tsv \
/path/to/A1_Inf_802/abundance.tsv \
/path/to/A1_MixFer/abundance.tsv \
/path/to/A1_MixInf/abundance.tsv \
/path/to/A2_Fer_801/abundance.tsv \
/path/to/A2_Fer_841/abundance.tsv \
/path/to/A2_Inf_001/abundance.tsv \
/path/to/A2_Inf_802/abundance.tsv \
/path/to/A2_MixFer/abundance.tsv \
/path/to/A2_MixInf/abundance.tsv \
/path/to/B1_Fer_801/abundance.tsv \
/path/to/B1_Fer_841/abundance.tsv \
/path/to/B1_Inf_001/abundance.tsv \
/path/to/B1_Inf_802/abundance.tsv \
/path/to/B1_MixFer/abundance.tsv \
/path/to/B1_MixInf/abundance.tsv \
/path/to/B2_Fer_801/abundance.tsv \
/path/to/B2_Fer_841/abundance.tsv \
/path/to/B2_Inf_001/abundance.tsv \
/path/to/B2_Inf_802/abundance.tsv \
/path/to/B2_MixFer/abundance.tsv \
/path/to/B2_MixInf/abundance.tsv \
/path/to/C1_Fer_801/abundance.tsv \
/path/to/C1_Fer_841/abundance.tsv \
/path/to/C1_Inf_001/abundance.tsv \
/path/to/C1_Inf_802/abundance.tsv \
/path/to/C1_MixFer/abundance.tsv \
/path/to/C1_MixInf/abundance.tsv \
/path/to/C2_Fer_801/abundance.tsv \
/path/to/C2_Fer_841/abundance.tsv \
/path/to/C2_Inf_001/abundance.tsv \
/path/to/C2_Inf_802/abundance.tsv \
/path/to/C2_MixFer/abundance.tsv \
/path/to/C2_MixInf/abundance.tsv \
/path/to/D_Fer_801/abundance.tsv \
/path/to/D_Fer_841/abundance.tsv \
/path/to/D_Inf_001/abundance.tsv \
/path/to/D_Inf_802/abundance.tsv \
/path/to/D_MixFer/abundance.tsv \
/path/to/D_MixInf/abundance.tsv \
/path/to/ALDH_Loser_1_minus/abundance.tsv \
/path/to/ALDH_Loser_1_plus/abundance.tsv \
/path/to/ALDH_Loser_7_minus/abundance.tsv \
/path/to/ALDH_Loser_7_plus/abundance.tsv \
/path/to/ALDH_Winner_2B_minus/abundance.tsv \
/path/to/ALDH_Winner_2B_plus/abundance.tsv \
/path/to/ALDH_Winner_6_minus/abundance.tsv \
/path/to/ALDH_Winner_6_plus/abundance.tsv \
/path/to/Ampullae_SB825FDA_Rep1/abundance.tsv \
/path/to/Ampullae_SB825FDA_Rep2/abundance.tsv \
/path/to/Ampullae_SB825FDA_Rep3/abundance.tsv \
/path/to/Ampullae_SB825FEA_Rep1/abundance.tsv \
/path/to/Ampullae_SB825FEA_Rep2/abundance.tsv \
/path/to/Ampullae_SB825FEA_Rep3/abundance.tsv \
/path/to/Ampullae_SB825FX_Rep1/abundance.tsv \
/path/to/Ampullae_SB825FX_Rep2/abundance.tsv \
/path/to/Ampullae_SB825FX_Rep3/abundance.tsv \
/path/to/Ampullae_SMGO_Rep2/abundance.tsv \
/path/to/Ampullae_SMGO_Rep3/abundance.tsv \
/path/to/Ampullae_SMGO_Technical_Rep1/abundance.tsv \
/path/to/Ampullae_SMGO_Technical_Rep2/abundance.tsv \
/path/to/Ampullae_SMH_Rep1/abundance.tsv \
/path/to/Ampullae_SMH_Rep2/abundance.tsv \
/path/to/Ampullae_SMH_Rep3/abundance.tsv \
/path/to/Blood_FL3/abundance.tsv \
/path/to/Blood_FLLC2/abundance.tsv \
/path/to/Blood_LFC4/abundance.tsv \
/path/to/IA6_79C_plus_rep1/abundance.tsv \
/path/to/IA6_79C_plus_rep2/abundance.tsv \
/path/to/IA6_79G_plus_rep1/abundance.tsv \
/path/to/IA6_79G_plus_rep2/abundance.tsv \
/path/to/IA6_81A_minu_rep1/abundance.tsv \
/path/to/IA6_81A_minu_rep2/abundance.tsv \
/path/to/IA6_81A_plus_rep1/abundance.tsv \
/path/to/IA6_81A_plus_rep2/abundance.tsv \
/path/to/IA6_81C_minu_rep1/abundance.tsv \
/path/to/IA6_81C_minu_rep2/abundance.tsv \
/path/to/IA6_81C_plus_rep1/abundance.tsv \
/path/to/IA6_81C_plus_rep2/abundance.tsv \
/path/to/Male_M001_A1/abundance.tsv \
/path/to/Male_M001_A2/abundance.tsv \
/path/to/Male_M001_B1/abundance.tsv \
/path/to/Male_M001_B2/abundance.tsv \
/path/to/Male_M001_C1/abundance.tsv \
/path/to/Male_M001_C2/abundance.tsv \
/path/to/Male_M001_D/abundance.tsv \
/path/to/Male_MADM_A1/abundance.tsv \
/path/to/Male_MADM_A2/abundance.tsv \
/path/to/Male_MADM_B1/abundance.tsv \
/path/to/Male_MADM_B2/abundance.tsv \
/path/to/Male_MADM_C1/abundance.tsv \
/path/to/Male_MADM_C2/abundance.tsv \
/path/to/Male_MADM_D/abundance.tsv \
/path/to/Male_MDR_A1/abundance.tsv \
/path/to/Male_MDR_A2/abundance.tsv \
/path/to/Male_MDR_B1/abundance.tsv \
/path/to/Male_MDR_B2/abundance.tsv \
/path/to/Male_MDR_C1/abundance.tsv \
/path/to/Male_MDR_C2/abundance.tsv \
/path/to/Male_MDR_D/abundance.tsv \
/path/to/Oocyte_with_follicles/abundance.tsv \
/path/to/Oocyte_without_follicles/abundance.tsv \
/path/to/WB_825_rep1/abundance.tsv \
/path/to/WB_825_rep2/abundance.tsv \
/path/to/WB_825_rep3/abundance.tsv \
/path/to/WB_SKI_rep1/abundance.tsv \
/path/to/WB_SKI_rep2/abundance.tsv \
/path/to/WB_SKI_rep3/abundance.tsv \
/path/to/WB_SMGO_rep1/abundance.tsv \
/path/to/WB_SMGO_rep2/abundance.tsv \
/path/to/WB_SMGO_rep3/abundance.tsv \
/path/to/WB_SMH_802/abundance.tsv \
/path/to/WB_SMH_rep1/abundance.tsv \
/path/to/WB_SMH_rep2/abundance.tsv \
/path/to/BSA_vasculature_cells_04/abundance.tsv \
/path/to/BSA_vasculature_cells_05/abundance.tsv \
/path/to/BSA_vasculature_cells_06/abundance.tsv \
/path/to/BSA_vasculature_cells_07/abundance.tsv \
/path/to/BSA_vasculature_cells_10/abundance.tsv \
/path/to/T0_rep1/abundance.tsv \
/path/to/T0_rep2/abundance.tsv \
/path/to/T0_rep3/abundance.tsv \
/path/to/T6_rep1/abundance.tsv \
/path/to/T6_rep2/abundance.tsv \
/path/to/T6_rep3/abundance.tsv \
/path/to/T18_rep1/abundance.tsv \
/path/to/T18_rep2/abundance.tsv \
/path/to/T18_rep3/abundance.tsv \
/path/to/T24_rep1/abundance.tsv \
/path/to/T24_rep2/abundance.tsv \
/path/to/T24_rep3/abundance.tsv


# Deactivate Conda environment
conda deactivate