# Load necessary libraries
library(ggplot2)
library(reshape2)
library(argparse)
library(readr)
library(purrr)

parser <- ArgumentParser()
parser$add_argument("--input", type = "character", nargs="+", help = "Input Spikein Count CSV file")
parser$add_argument("--output", type = "character", help = "Output PNG file")

# Parsing arguments
args <- parser$parse_args()

validate_data <- function(counts_data) {
  if (nrow(counts_data) == 0) {
    stop("No spikein count data was detected.")
  }
}

melt_data <- function(counts_data) {

  # Transform the data to logaritmic
  counts_data$Count <- log10(counts_data$Count + 1)

  # Transform the data from long to wide format
  wide_data <- dcast(counts_data, Sample ~ SpikeinID, value.var = "Count")

  # Melt the data for ggplot2
  melted_data <- melt(wide_data, id.vars = "Sample")

  # Return melted data
  return (melted_data)
}

# Plot where spikeins (x-axis) are detected by sample (y-axis)
plot_spikein_detection_heatmap_by_sampleid <- function(melted_data) {

  # Generate the heatmap with a white background
  g <- ggplot(melted_data, aes(x = variable, y = Sample, fill = value)) +
    geom_tile() +
    scale_fill_gradient(name = "Log10(Count + 1)", low = "blue", high = "red") + # Indicate log scale in legend
    labs(x = "Spike-in ID", y = "Sample", fill = "Count") +
    theme_minimal(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))

  return (g)
}

# Load the data
counts_data <- args$input |> map_dfr(read.csv)
print(counts_data)
validate_data(counts_data)
melted_data <- melt_data(counts_data)

# Save the heatmap with a white background
g <- plot_spikein_detection_heatmap_by_sampleid(melted_data)
ggsave(args$output, plot = g, bg = "white", width = 10, height = 8, dpi = 300)