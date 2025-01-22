#!/usr/bin/env Rscript

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_barplot.R <input_csv> <output_png> <ratio>")
}

input_csv <- args[1]
output_png <- args[2]
ratio <- as.numeric(args[3])

# Read the CSV data
data <- read.csv(input_csv)

# Filter out rows where the proportion is zero
data <- data[data$Proportion > 0, ]

# Add a column to classify motifs for coloring
data$Type <- ifelse(
  data$Motif == "LowComplexity", "LowComplexity",
  ifelse(nchar(as.character(data$Motif)) == 2, "DiNucleotide", "TriNucleotide")
)

# Define custom colors for the motif types
colors <- c("LowComplexity" = "red", "DiNucleotide" = "blue", "TriNucleotide" = "lightblue")

# Create the barplot
library(ggplot2)
p <- ggplot(data, aes(x = reorder(Motif, -Proportion), y = Proportion, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(
    title = sprintf("Proportion of reads with at least %.1f%% of given motif", ratio),
    x = "Motif",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  ) +
  coord_cartesian(ylim = c(0, 0.05))  # Keep bars beyond 0.05 but limit the y-axis display

# Add labels for bars exceeding 0.05
p <- p + geom_text(
  aes(
    label = ifelse(Proportion > 0.05, sprintf("%.4f", Proportion), ""),
    y = ifelse(Proportion > 0.05, 0.05, Proportion)
  ),
  vjust = 0.5,
  color = "black",
  size = 4
)

# Save the plot to a PNG file with a fully white background
ggsave(output_png, plot = p, width = 12, height = 8, bg = "white")
