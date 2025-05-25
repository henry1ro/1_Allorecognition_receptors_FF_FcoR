#Barplot showing the number of Fester and Fester coreceptor genes per individual

# Load necessary libraries
library(ggplot2)
library(tidyr)

# Input the new mean data from the image
expression_data <- data.frame(
  Individual = c('F841','F801','Infer001','Infer802','M001','MADM','MDR','Ampullae-SB825FDA','Ampullae-SB825FEA','Ampullae-SB825FX','Ampullae-SMGO','Ampullae-SMH','WB-825','WB-SKI','WB-SMH','WB-SMGO','Aging-SB825','Aging-SMH','Aging-SMGO','WBR'),
  Fester = c(9,11,10,10,10,11,9,9,10,11,7,10,10,10,11,7,10,13,7,7), #number of fester genes per individual
  Fester_coreceptor = c(9,11,10,12,13,13,12,7,8,6,4,10,28,19,26,18,15,17,20,15) #number of fester coreceptor genes per individual
)

# Reshape data to long format
expression_data_long <- gather(expression_data, key = "Gene", value = "RNA_Expression", -Individual)

# Reorder levels of the "Gene" variable
expression_data_long$Gene <- factor(expression_data_long$Gene, levels = c("Fester", "Fester_coreceptor"))

# Convert "Individual" to factor with desired levels
expression_data_long$Individual <- factor(expression_data_long$Individual, levels = c('F841','F801','Infer001','Infer802','M001','MADM','MDR',
                 'Ampullae-SB825FDA','Ampullae-SB825FEA','Ampullae-SB825FX','Ampullae-SMGO','Ampullae-SMH',
                 'WB-825','WB-SKI','WB-SMH','WB-SMGO','Aging-SB825','Aging-SMH','Aging-SMGO','WBR',
                 'MixFer','MixInfer','ALDH','BSA','IA6','Blood')) #individuals

# Create the histogram with larger x-axis label size
histogram <- ggplot(expression_data_long, aes(x = Individual, y = RNA_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Number of Fester and Fester-coreceptor genes per individual",
       x = "Individuals",
       y = "Number of genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),  # Increased size
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.key.size = unit(0.6, "cm"),
        legend.key = element_rect(color = "black", size = 0.2),
        legend.position = c(0.88, 0.85)) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30))  # Adjust y-axis breaks and limits

# Print the histogram
print(histogram)
