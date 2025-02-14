#!/usr/bin/env Rscript

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: generate_barplot.R <input_csv> <output_png> <ratio>")
}

input_csv <- args[1]
output_png <- args[2]
ratio <- as.numeric(args[3])

# Read the CSV data
data_f <- read.csv(input_csv)

# Filter out rows where the proportion is zero
data <- data_f[data_f$Proportion > 1.0, ]

# Add a column to classify motifs for coloring
data$Type <- ifelse(
  data$Motif == "LowComplexity", "LowComplexity",
  ifelse(nchar(as.character(data$Motif)) == 2, "DiNucleotide", "TriNucleotide")
)

# Assign dynamic colors to "LowComplexity" based on its value
data$LowComplexityColor <- ifelse(
  data$Motif == "LowComplexity" & data$Proportion > 15, "red",
  ifelse(data$Motif == "LowComplexity" & data$Proportion > 5, "orange", 
         ifelse(data$Motif == "LowComplexity", "green", NA))
)

# Extract the LowComplexity proportion value
low_complexity_value <- data_f$Proportion[data_f$Motif == "LowComplexity"]

# Assign colors for other types
data$TypeColor <- ifelse(
  data$Type == "DiNucleotide" & data$Proportion < low_complexity_value & data$Proportion > 5 & !(data$Motif %in% c("AA", "TT", "CC", "GG")), 
  "red", 
  ifelse(data$Type == "DiNucleotide", "blue", "lightblue")
)

# Combine colors into a single column for plotting
data$FinalColor <- ifelse(!is.na(data$LowComplexityColor), data$LowComplexityColor, data$TypeColor)

# Create the barplot
library(ggplot2)
p <- ggplot(data, aes(x = reorder(Motif, -Proportion), y = Proportion, fill = FinalColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() + 
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
  coord_cartesian(ylim = c(0, 50))

# Add labels for bars exceeding 50
p <- p + geom_text(
  aes(
    label = sprintf("%.1f", Proportion),
    y = ifelse(Proportion > 50, 50, Proportion)
  ),
  vjust = -0.1,
  color = "black",
  size = 4
)

# Save the plot to a PNG file with a fully white background
ggsave(output_png, plot = p, width = 12, height = 8, bg = "white")
