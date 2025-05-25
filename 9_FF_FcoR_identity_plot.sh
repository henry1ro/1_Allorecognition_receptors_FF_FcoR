#Generate protein identity plots for Fester and Fester coreceptor proteins

#Calculate an protein identity matrix using software such as BioEdit

# Load necessary packages
library(ggplot2)
library(reshape2)
library(readr)
library(viridis)

# Read the data from the text file
identity_matrix <- read_delim("/path/to/matrix_fester.txt", delim = "\t", col_names = FALSE) #Identity matrix for Fester and Fester coreceptor proteins

# Convert to data frame
identity_matrix <- as.data.frame(identity_matrix)

# Set the first row as column names
colnames(identity_matrix) <- identity_matrix[1,]
identity_matrix <- identity_matrix[-1,]

# Set the first column as row names
rownames(identity_matrix) <- identity_matrix[[1]]
identity_matrix <- identity_matrix[-1]

# Convert data to numeric
identity_matrix <- as.data.frame(lapply(identity_matrix, as.numeric))

# Melt the matrix into long format
identity_matrix_long <- melt(as.matrix(identity_matrix), varnames = c("Protein1", "Protein2"), value.name = "Identity")

# Convert Protein1 and Protein2 to factors to preserve categorical ordering
identity_matrix_long$Protein1 <- factor(identity_matrix_long$Protein1, levels = rownames(identity_matrix))
identity_matrix_long$Protein2 <- factor(identity_matrix_long$Protein2, levels = colnames(identity_matrix))

# Create the heatmap with improved aesthetics
heatmap_plot <- ggplot(identity_matrix_long, aes(x = Protein1, y = Protein2, fill = Identity)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", name = "Identity (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, margin = margin(t = 5)),
    axis.text.y = element_text(size = 8, margin = margin(r = 5)),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "Fester Protein Identity", #Title
    x = "Protein",
    y = "Protein"
  )

# Display the plot
print(heatmap_plot)

# Save the plot to a file
ggsave("protein_identity_heatmap.png", plot = heatmap_plot, width = 10, height = 10, dpi = 300)



#Boxplot comparing the identities of FF genes vs FF1 alleles (or FcoR genes vs FcoR7 alelles)

# Install and load required packages
install.packages("ggplot2")
install.packages("reshape2")
library(ggplot2)
library(reshape2)

# Read the data from the CSV file
data <- read.csv("/path/to/data_fester.csv") #File with two columns: the first column contains the identities for genes, and the second column contains the identities for alleles

# Melt the data for ggplot2
data_melted <- melt(data)

# Replace the labels in the melted data
data_melted$variable <- gsub("\\.", " ", data_melted$variable)  # Remove dots in variable names

# Create the boxplot with y-axis divisions of 0.1
ggplot(data_melted, aes(x=variable, y=value)) +
  geom_boxplot(width=0.7,fill=c("blue", "red")) +
  ylab("Identity") +
  xlab("") +
  scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
  theme_classic() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(face = "bold", color = "black", size = 14))  # Adjusted x-axis text

