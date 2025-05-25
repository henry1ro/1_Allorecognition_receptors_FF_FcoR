#Heatmap showing the presence or absence of Fester and Fester coreceptor genes from transcriptome screening

library(ComplexHeatmap)

# Load data
data <- read.csv("/path/to/FF_FcoR_transcriptomic_data.csv", row.names = 1) #matrix indicating the presence (1) and absence (0) of FF and FcoR genes

# Convert data frame to numeric matrix
data <- as.matrix(data)

# Transpose the data to swap x and y axes
data_transposed <- t(data)

# Create heatmap with row dendrogram on the right and transposed axes
Heatmap(
  data_transposed,
  cluster_rows = TRUE,         # Cluster rows (which are now the original columns)
  cluster_columns = TRUE,      # Cluster columns (which are now the original rows)
  show_row_dend = TRUE,        # Show row dendrogram
  show_column_dend = FALSE,    # Hide column dendrogram
  row_dend_side = "right",     # Position row dendrogram to the right
  column_title = "Presence of FF and FcoR in transcriptomes", # Title of the heatmap
  col = colorRampPalette(c("white", "black"))(100), # Black and white color scheme
  show_heatmap_legend = FALSE, # Remove the heatmap legend
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "black", fill = fill, lwd = 0.5))  # Black borders inside cells
  },
  row_names_gp = gpar(fontsize = 8),  # Change the font size of row names (genes)
  column_names_gp = gpar(fontsize = 8),  # Change the font size of column names
  column_title_gp = gpar(fontsize = 10)  # Change the font size of the column title
)
