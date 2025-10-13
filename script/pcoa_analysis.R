#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(svglite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript pcoa_analysis.R <input_csv> <output_svg>")
}

input_file <- args[1]
output_file <- args[2]

# Read abundance table
otu_data <- read.csv(input_file, check.names = FALSE)

# Check if table has taxlevel, taxid, taxa columns (Kraken format)
if (ncol(otu_data) > 3 && all(c("taxlevel", "taxid", "taxa") %in% colnames(otu_data))) {
  # Extract taxa names and make unique if needed
  taxa_names <- otu_data$taxa
  if (any(duplicated(taxa_names))) {
    taxa_names <- make.unique(taxa_names, sep = "_")
  }

  # Keep only sample columns (skip taxlevel, taxid, taxa)
  sample_cols <- setdiff(colnames(otu_data), c("taxlevel", "taxid", "taxa"))
  otu_data <- otu_data[, sample_cols, drop = FALSE]
  rownames(otu_data) <- taxa_names
} else {
  # Assume first column is taxa names
  taxa_names <- otu_data[[1]]
  if (any(duplicated(taxa_names))) {
    taxa_names <- make.unique(taxa_names, sep = "_")
  }
  rownames(otu_data) <- taxa_names
  otu_data <- otu_data[, -1, drop = FALSE]
}

# Convert to numeric and handle missing values
otu_data[] <- lapply(otu_data, function(x) as.numeric(as.character(x)))

# Transpose data (samples as rows, taxa as columns)
otu_data <- data.frame(t(otu_data))
otu_data[is.na(otu_data)] <- 0

# Clean up sample names - remove kraken/bracken suffixes
rownames(otu_data) <- gsub("_kraken2_bracken_.*\\.report$", "", rownames(otu_data))

# Simplify sample names: W127-25PE0.6-2-17 -> PE0.6-2-17 (remove prefix like Q31-25 or W127-30)
rownames(otu_data) <- gsub("^[A-Z]\\d+-\\d+", "", rownames(otu_data))

# Calculate Bray-Curtis distance
bray_dist <- vegdist(otu_data, method = 'bray')

# Perform PCoA
pcoa_result <- cmdscale(bray_dist, k = min(nrow(otu_data) - 1, 5), eig = TRUE)

# Calculate variance explained by each axis
variance_exp <- pcoa_result$eig / sum(pcoa_result$eig)

# Prepare data for plotting
site_scores <- data.frame(pcoa_result$point[, 1:2])
colnames(site_scores) <- c("PCoA1", "PCoA2")
site_scores$Sample <- rownames(site_scores)

# Extract group information (remove replicate number -1-, -2-, -3-)
# Pattern: PE0.6-2-17 -> PE0.6-17 (remove replicate number for grouping)
site_scores$Group <- gsub("^(.+)-[123]-(\\d+)$", "\\1-\\2", site_scores$Sample)

# Create color palette for groups
n_groups <- length(unique(site_scores$Group))
group_colors <- scales::hue_pal()(n_groups)
names(group_colors) <- sort(unique(site_scores$Group))

# Create PCoA plot
pcoa_plot <- ggplot(data = site_scores, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Group), size = 4, shape = 16) +
  scale_color_manual(values = group_colors) +
  geom_text(aes(label = Sample), vjust = -1, size = 2.5, color = "black") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = 'gray', linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'white'),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, color = 'gray60', linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, color = 'gray60', linewidth = 0.5, linetype = "dashed") +
  labs(
    x = paste0("PCoA1 (", round(variance_exp[1] * 100, 2), "%)"),
    y = paste0("PCoA2 (", round(variance_exp[2] * 100, 2), "%)")
  )

# Save as SVG (editable)
ggsave(output_file, pcoa_plot, width = 8, height = 6, device = "svg")

cat("PCoA plot saved to:", output_file, "\n")
