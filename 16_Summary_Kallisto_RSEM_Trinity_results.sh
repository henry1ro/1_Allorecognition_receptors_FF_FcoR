#Summary of the presence/absence of FF/FcoR gene pairs across Trinity, Kallisto, and RSEM methods

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ggtext)

# Read the tab-delimited .txt file
data <- read.table("/path/to/table.txt", sep = "\t", header = TRUE) #File indicating the presence (1) and absence (0) of FF/FcoR gene pairs across Trinity, Kallisto, and RSEM methods

# Combine the first two columns to create unique row identifiers
data$Sample <- paste(data$Fester, data$FcoR, sep = "_")

# Remove the original first two columns from the dataset for melting
data_for_melt <- data[, !(colnames(data) %in% c("Fester", "FcoR"))]

# Melt the data for ggplot2
data_melt <- melt(data_for_melt, id.vars = "Sample", variable.name = "Method", value.name = "Presence")

# Format the sample names: highlight FF in blue and FcoR in red
data_melt$Sample <- gsub(
  "(FF\\d+)", 
  "<span style='color:blue;'>\\1</span>", 
  data_melt$Sample
)
data_melt$Sample <- gsub(
  "(FcoR\\d+)", 
  "<span style='color:red;'>\\1</span>", 
  data_melt$Sample
)

# Ensure the sample order is preserved
data_melt$Sample <- factor(data_melt$Sample, levels = unique(data_melt$Sample))

# Create the heatmap plot
ggplot(data_melt, aes(x = Method, y = Sample, fill = as.factor(Presence))) +
  geom_tile(color = "black") +
  scale_fill_manual(
    values = c("0" = "white", "1" = "black"),
    guide = "none"  # Remove the legend
  ) +
  labs(x = "Methods", y = "Gene pairs", title = "FF and FcoR gene pairs") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # Center and increase title size
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Center and increase x-axis text size
    axis.text.y = element_markdown(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14)   # Increase y-axis title size
  )
