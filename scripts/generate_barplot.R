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

# Assign dynamic colors to "LowComplexity" based on its value
data$LowComplexityColor <- ifelse(
  data$Motif == "LowComplexity" & data$Proportion > 10, "red",
  ifelse(data$Motif == "LowComplexity" & data$Proportion > 5, "orange", 
         ifelse(data$Motif == "LowComplexity", "green", NA))
)

# Assign colors for other types
data$TypeColor <- ifelse(data$Type == "DiNucleotide", "blue", "lightblue")

# Combine colors into a single column for plotting
data$FinalColor <- ifelse(!is.na(data$LowComplexityColor), data$LowComplexityColor, data$TypeColor)

# Create the barplot
library(ggplot2)
p <- ggplot(data, aes(x = reorder(Motif, -Proportion), y = Proportion, fill = FinalColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  # Use the exact colors specified in the FinalColor column
  labs(
    title = sprintf("Proportion of reads with at least %.0f%% of the specified motif", ratio),
    x = "Motif",
    y = "Proportion (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  ) +
  coord_cartesian(ylim = c(0, 0.5))  # Keep bars beyond 0.05 but limit the y-axis display

# Add labels for bars exceeding 0.5
p <- p + geom_text(
  aes(
    label = ifelse(Proportion > 0.1, 
               sprintf("%.2f", Proportion), 
               ifelse(Proportion > 0.05, 
                      sprintf("%.3f", Proportion), 
                      "")),
    y = ifelse(Proportion > 0.5, 0.5, Proportion)
  ),
  vjust = -0.1,
  color = "black",
  size = 4
)

# Save the plot to a PNG file with a fully white background
ggsave(output_png, plot = p, width = 12, height = 8, bg = "white")
