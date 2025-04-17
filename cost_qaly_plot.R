# Script to create a scatterplot of costs vs QALYs for ESD and Surgery

# Load required packages
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require(reshape2)) {
  install.packages("reshape2")
  library(reshape2)
}

# Read in the simulation results
results_dir <- "./results"
results_file <- file.path(results_dir, "simulation_results_live.csv")

if (!file.exists(results_file)) {
  stop("Results file not found! Please check the path:", results_file)
}

# Read the data
sim_results <- read.csv(results_file)

# Create a long-format data frame for the scatterplot
# We'll create two datasets - one for ESD and one for Surgery
esd_data <- data.frame(
  Strategy = rep("ESD", nrow(sim_results)),
  QALYs = sim_results$esd_qalys,
  Cost = sim_results$esd_cost,
  Simulation = sim_results$simulation
)

surgery_data <- data.frame(
  Strategy = rep("Surgery", nrow(sim_results)),
  QALYs = sim_results$surgery_qalys,
  Cost = sim_results$surgery_cost,
  Simulation = sim_results$simulation
)

# Combine the datasets
plot_data <- rbind(esd_data, surgery_data)

# Create the scatterplot
p <- ggplot(plot_data, aes(x = QALYs, y = Cost, color = Strategy)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("ESD" = "blue", "Surgery" = "red")) +
  labs(
    title = "Cost vs QALYs by Treatment Strategy",
    x = "Quality-Adjusted Life Years (QALYs)",
    y = "Cost ($)",
    color = "Strategy"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Show the plot in the R console
print(p)

# Save the plot to a file
ggsave(file.path(results_dir, "cost_qaly_scatterplot.png"), p, width = 10, height = 8)

# Create an additional plot with convex hulls around each strategy group
p_hull <- p + 
  geom_polygon(
    data = plot_data,
    aes(x = QALYs, y = Cost, fill = Strategy),
    alpha = 0.2,
    stat = "chull"
  ) +
  scale_fill_manual(values = c("ESD" = "blue", "Surgery" = "red"))

# Save the second plot to a file
ggsave(file.path(results_dir, "cost_qaly_scatterplot_with_hulls.png"), p_hull, width = 10, height = 8)

# Create a boxplot for QALYs by strategy
p_boxplot <- ggplot(plot_data, aes(x = Strategy, y = QALYs, fill = Strategy)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ESD" = "blue", "Surgery" = "red")) +
  labs(
    title = "Distribution of QALYs by Treatment Strategy",
    x = "",
    y = "Quality-Adjusted Life Years (QALYs)",
    fill = "Strategy"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Save the boxplot to a file
ggsave(file.path(results_dir, "qalys_boxplot_by_strategy.png"), p_boxplot, width = 8, height = 6)

# Create a boxplot for Costs by strategy
p_cost_boxplot <- ggplot(plot_data, aes(x = Strategy, y = Cost, fill = Strategy)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ESD" = "blue", "Surgery" = "red")) +
  labs(
    title = "Distribution of Costs by Treatment Strategy",
    x = "",
    y = "Cost ($)",
    fill = "Strategy"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Save the cost boxplot to a file
ggsave(file.path(results_dir, "cost_boxplot_by_strategy.png"), p_cost_boxplot, width = 8, height = 6)

cat("Plots have been saved to the", results_dir, "directory.\n")
cat("Files created:\n")
cat("1. cost_qaly_scatterplot.png\n")
cat("2. cost_qaly_scatterplot_with_hulls.png\n")
cat("3. qalys_boxplot_by_strategy.png\n")
cat("4. cost_boxplot_by_strategy.png\n")
