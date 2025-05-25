#Mapping raw reads from vasculature transcriptomes to allorecognition genes

#!/bin/bash -l
# ask for 30 cores on three nodes
#SBATCH --nodes=3 --ntasks-per-node=30
#SBATCH --time=700:00:00

# Activate Conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate trin

align_and_estimate_abundance.pl \
--est_method kallisto \
--samples_file /path/to/SampleFile_TotalBotryllus_7_vas.txt \  #sample file that specifies the path to the vasculature raw reads
--transcripts /path/to/6_NK-cell-genes.fasta \   #allorecognition or ITAM/ITIM genes
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
--est_method kallisto \
--out_prefix all \
--gene_trans_map none \
--name_sample_by_basedir \
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
/path/to/Ampullae_SMH_Rep3/abundance.tsv

# Deactivate Conda environment
conda deactivate

#Barplot for allorecognition ligands and ITAM/ITIM genes

# Load necessary libraries
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(writexl)

# --- Data Loading and Initial Processing ---
# Read the Excel file
data_nk <- read_excel("/path/to/nk_vas.xlsx") #path to mapping results of allorecognition genes

# Rename the first column to Gene
colnames(data_nk)[1] <- "Gene"

# Reshape the data
long_data_nk <- data_nk %>%
  pivot_longer(cols = -Gene,
               names_to = c("Individual", "Replicate"),
               names_pattern = "(.*)_(Rep[1-4])",
               values_to = "Expression")

# --- ANOVA Analysis ---
# Run ANOVA for each gene
anova_results_nk <- long_data_nk %>%
  group_by(Gene) %>%
  group_modify(~ {
    result <- aov(Expression ~ Individual, data = .x)
    tidy_result <- broom::tidy(result)
    tidy_result
  }) %>%
  filter(term == "Individual") %>%
  mutate(adjusted_p_value = p.adjust(p.value, method = "BH"))

# Handle NaN values in adjusted_p_value
anova_results_nk <- anova_results_nk %>%
  mutate(adjusted_p_value = ifelse(is.nan(adjusted_p_value), 1, adjusted_p_value))

# Create significance annotations
anova_results_summary_nk <- anova_results_nk %>%
  mutate(Significance = case_when(
    adjusted_p_value < 0.001 ~ "***",
    adjusted_p_value < 0.01  ~ "**",
    adjusted_p_value < 0.05  ~ "*",
    TRUE                    ~ ""
  )) %>%
  select(Gene, term, df, sumsq, meansq, statistic, p.value, adjusted_p_value, Significance)

# --- Summary Statistics for Plotting ---
# Calculate mean expression and standard deviation for each gene and individual
mean_expression_nk <- long_data_nk %>%
  group_by(Gene, Individual) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    SD_Expression = sd(Expression, na.rm = TRUE)
  ) %>%
  mutate(SD_Expression = pmax(SD_Expression, 0)) # Ensure SD is not below zero

# Filter to only include rows where SD > 0
mean_expression_filtered_nk <- mean_expression_nk %>%
  filter(SD_Expression > 0)

# Calculate the maximum expression level for each gene
max_expression_nk <- mean_expression_nk %>%
  group_by(Gene) %>%
  summarise(Max_Expression = max(Mean_Expression, na.rm = TRUE))

# Merge data for plotting
plot_data_nk <- mean_expression_filtered_nk %>%
  left_join(anova_results_summary_nk, by = "Gene") %>%
  left_join(max_expression_nk, by = "Gene")

# --- Define Gene Order ---
desired_gene_order <- c("BHF", "FuHc-Sec", "FuHc-Tm", "HSP40L", "SYK", "SHP2", "SHIP", "PI3K", "Src", "Grb2", "VAV", "Sos", "Shc", "PLC-g", "PKC", "NFAT", "NF-kB", "Ras", "Raf", "Rac", "PAK1", "MEK", "ERK", "CnA", "CnB", "E3-ligase")

# Filter plot data to include only the desired genes and set the factor level
plot_data_nk_filtered <- plot_data_nk %>%
  filter(Gene %in% desired_gene_order) %>%
  mutate(Gene = factor(Gene, levels = desired_gene_order)) %>%
  group_by(Gene) %>%
  mutate(Max_Expression_Annotated = max(Mean_Expression + SD_Expression, na.rm = TRUE) + 0.1)

# --- Define Color Palette ---
custom_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a")

# --- Create and Print Bar Plot with Significance ---
plot_nk <- ggplot(plot_data_nk_filtered, aes(x = Gene, y = Mean_Expression, fill = Individual)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = pmax(Mean_Expression - SD_Expression, 0),
                    ymax = Mean_Expression + SD_Expression),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5,
                color = "black",
                linetype = "solid") +
  geom_text(aes(x = Gene, y = Max_Expression_Annotated, label = Significance),
            vjust = -0.5, size = 6, color = "black") +
  scale_fill_manual(values = custom_palette) +
  theme_classic() +
  labs(title = "Expression of ITAM/ITIM signal transduction genes in Botryllus vasculature",
       x = "Genes",
       y = "RNA expression (TMM value)",
       fill = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  coord_cartesian(clip = 'off')

# Print the plot
print(plot_nk)

# --- Export ANOVA Results for the desired genes ---
anova_results_summary_nk_filtered <- anova_results_summary_nk %>%
  filter(Gene %in% desired_gene_order)

write_xlsx(anova_results_summary_nk_filtered, "NK_anova.xlsx")


#Barplots for FF and FcoR genes

# Load necessary libraries
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(writexl)

# --- Data Loading and Initial Processing ---
# Read the Excel files
data_FF <- read_excel("/path/to/fester.xlsx") #path to mapping results of FF genes
data_FcoR <- read_excel("/path/to/FcoR.xlsx") #path to mapping results of FcoR genes

# Rename the first column to Gene for both datasets
colnames(data_FF)[1] <- "Gene"
colnames(data_FcoR)[1] <- "Gene"

# Reshape the data for both datasets
long_data_FF <- data_FF %>%
  pivot_longer(cols = -Gene,
               names_to = c("Individual", "Replicate"),
               names_pattern = "(.*)_(Rep[1-4])",
               values_to = "Expression")

long_data_FcoR <- data_FcoR %>%
  pivot_longer(cols = -Gene,
               names_to = c("Individual", "Replicate"),
               names_pattern = "(.*)_(Rep[1-4])",
               values_to = "Expression")

# Combine the long data
long_data <- rbind(long_data_FF, long_data_FcoR)

# --- ANOVA Analysis ---
# Run ANOVA for each gene
anova_results <- long_data %>%
  group_by(Gene) %>%
  group_modify(~ {
    result <- aov(Expression ~ Individual, data = .x)
    tidy_result <- broom::tidy(result)
    tidy_result
  }) %>%
  filter(term == "Individual") %>%
  mutate(adjusted_p_value = p.adjust(p.value, method = "BH"))

# Handle NaN values in adjusted_p_value
anova_results <- anova_results %>%
  mutate(adjusted_p_value = ifelse(is.nan(adjusted_p_value), 1, adjusted_p_value))

# Create significance annotations
anova_results_summary <- anova_results %>%
  mutate(Significance = case_when(
    adjusted_p_value < 0.001 ~ "***",
    adjusted_p_value < 0.01  ~ "**",
    adjusted_p_value < 0.05  ~ "*",
    TRUE                    ~ ""
  )) %>%
  select(Gene, term, df, sumsq, meansq, statistic, p.value, adjusted_p_value, Significance)

# --- Summary Statistics for Plotting ---
# Calculate mean expression and standard deviation for each gene and individual
mean_expression <- long_data %>%
  group_by(Gene, Individual) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    SD_Expression = sd(Expression, na.rm = TRUE)
  ) %>%
  mutate(SD_Expression = pmax(SD_Expression, 0)) # Ensure SD is not below zero

# Filter to only include rows where SD > 0
mean_expression_filtered <- mean_expression %>%
  filter(SD_Expression > 0)

# Calculate the maximum expression level for each gene
max_expression <- mean_expression %>%
  group_by(Gene) %>%
  summarise(Max_Expression = max(Mean_Expression, na.rm = TRUE))

# Merge data for plotting
plot_data <- mean_expression_filtered %>%
  left_join(anova_results_summary, by = "Gene") %>%
  left_join(max_expression, by = "Gene")

# --- Define Gene Orders ---
ff_gene_order <- c("FF1", "FF3", "FF4", "FF7", "FF9", "FF10", "FF16", "FF24", "FF25", "FF28", "FF12", "FF13", "FF14", "FF15", "FF17", "FF19", "FF20", "FF44", "FF2", "FF5", "FF6", "FF8", "FF11", "FF23", "FF26", "FF29", "FF31", "FF35", "FF37", "FF40", "FF42")
fcor_gene_order_all <- c("FcoR7", "FcoR1", "FcoR12", "FcoR23", "FcoR22", "FcoR18", "FcoR46", "FcoR6", "FcoR35", "FcoR21", "FcoR3", "FcoR4", "FcoR16", "FcoR47", "FcoR63", "FcoR64", "FcoR5", "FcoR10", "FcoR17", "FcoR33", "FcoR34", "FcoR38", "FcoR39", "FcoR40", "FcoR43", "FcoR48", "FcoR2", "FcoR8", "FcoR9", "FcoR19", "FcoR20", "FcoR24", "FcoR26", "FcoR28", "FcoR31", "FcoR36", "FcoR44", "FcoR45", "FcoR50", "FcoR51", "FcoR52", "FcoR53", "FcoR54", "FcoR55", "FcoR57", "FcoR58", "FcoR59", "FcoR60", "FcoR61", "FcoR62", "FcoR65", "FcoR69", "FcoR70", "FcoR73", "FcoR75")

# Exclude FcoR50 from the FcoR gene order
fcor_gene_order <- fcor_gene_order_all[!fcor_gene_order_all %in% "FcoR50"]

# Filter and order FF genes
plot_data_ff <- plot_data %>%
  filter(grepl("^FF", Gene)) %>%
  mutate(Gene = factor(Gene, levels = ff_gene_order)) %>%
  group_by(Gene) %>%
  mutate(Max_Expression_Annotated = max(Mean_Expression + SD_Expression, na.rm = TRUE) + 0.1)

# Filter and order FcoR genes, excluding FcoR50
plot_data_fcor <- plot_data %>%
  filter(grepl("^FcoR", Gene) & Gene != "FcoR50") %>%
  mutate(Gene = factor(Gene, levels = fcor_gene_order)) %>%
  group_by(Gene) %>%
  mutate(Max_Expression_Annotated = max(Mean_Expression + SD_Expression, na.rm = TRUE) + 0.1)

# --- Define Color Palette ---
custom_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a")

# --- Create and Print Bar Plots with Significance ---
plot_ff <- ggplot(plot_data_ff, aes(x = Gene, y = Mean_Expression, fill = Individual)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = pmax(Mean_Expression - SD_Expression, 0),
                    ymax = Mean_Expression + SD_Expression),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5,
                color = "black",
                linetype = "solid") +
  geom_text(aes(x = Gene, y = Max_Expression_Annotated, label = Significance),
            vjust = -0.5, size = 6, color = "black") + # Changed color to "black"
  scale_fill_manual(values = custom_palette) +
  theme_classic() +
  labs(title = "Expression of FF genes in Botryllus vasculature",
       x = "Genes",
       y = "RNA expression (TMM value)",
       fill = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  coord_cartesian(clip = 'off')

plot_fcor <- ggplot(plot_data_fcor, aes(x = Gene, y = Mean_Expression, fill = Individual)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = pmax(Mean_Expression - SD_Expression, 0),
                    ymax = Mean_Expression + SD_Expression),
                position = position_dodge(width = 0.8),
                width = 0.2,
                size = 0.5,
                color = "black",
                linetype = "solid") +
  geom_text(aes(x = Gene, y = Max_Expression_Annotated, label = Significance),
            vjust = -0.5, size = 6, color = "black") + # Changed color to "black"
  scale_fill_manual(values = custom_palette) +
  theme_classic() +
  labs(title = "Expression of FcoR genes in Botryllus vasculature",
       x = "Genes",
       y = "RNA expression (TMM value)",
       fill = "Individual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  coord_cartesian(clip = 'off')

# Print the plots
print(plot_ff)
print(plot_fcor)

# --- Export ANOVA Results, excluding FcoR50 ---
anova_results_summary_filtered <- anova_results_summary %>%
  filter(Gene != "FcoR50")

write_xlsx(anova_results_summary_filtered, "FF_FcoR_anova.xlsx")






