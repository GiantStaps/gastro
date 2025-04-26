library(jsonlite)
library(tidyverse)  # for dplyr, tidyr, purrr
library(ggplot2)

# 1. Load and transform the JSON data
process_json_results <- function(json_path) {
  # Read the JSON file
  results <- fromJSON(json_path)
  
  # Create an empty list to store the data
  data_list <- list()
  
  # Loop through replications
  for (rep_id in names(results)) {
    rep_data <- results[[rep_id]]
    
    # Loop through cycles
    for (cycle_id in names(rep_data)) {
      cycle_data <- rep_data[[cycle_id]]
      
      # Process ESD data
      if (!is.null(cycle_data$esd)) {
        esd_states <- cycle_data$esd$state_counts
        
        # Convert state counts to long format
        for (state in names(esd_states)) {
          data_list[[length(data_list) + 1]] <- list(
            replication = as.integer(rep_id),
            cycle = as.integer(cycle_id),
            strategy = "esd",
            state = state,
            count = esd_states[[state]],
            total_cost = cycle_data$esd$cost$total$avg,
            total_cost_sd = cycle_data$esd$cost$total$std,
            surveillance_cost = cycle_data$esd$cost$surveillance$avg,
            surveillance_cost_sd = cycle_data$esd$cost$surveillance$std,
            qaly = cycle_data$esd$qaly$avg,
            qaly_sd = cycle_data$esd$qaly$std
          )
        }
      }
      
      # Process Surgery data
      if (!is.null(cycle_data$surgery)) {
        surgery_states <- cycle_data$surgery$state_counts
        
        # Convert state counts to long format
        for (state in names(surgery_states)) {
          data_list[[length(data_list) + 1]] <- list(
            replication = as.integer(rep_id),
            cycle = as.integer(cycle_id),
            strategy = "surgery",
            state = state,
            count = surgery_states[[state]],
            total_cost = cycle_data$surgery$cost$total$avg,
            total_cost_sd = cycle_data$surgery$cost$total$std,
            surveillance_cost = cycle_data$surgery$cost$surveillance$avg,
            surveillance_cost_sd = cycle_data$surgery$cost$surveillance$std,
            qaly = cycle_data$surgery$qaly$avg,
            qaly_sd = cycle_data$surgery$qaly$std
          )
        }
      }
      
      # Store ICER if available
      if (!is.null(cycle_data$icer) && !is.null(cycle_data$icer$avg)) {
        data_list[[length(data_list) + 1]] <- list(
          replication = as.integer(rep_id),
          cycle = as.integer(cycle_id),
          strategy = "icer",
          state = "ICER",
          count = NA,
          total_cost = NA,
          total_cost_sd = NA,
          surveillance_cost = NA,
          surveillance_cost_sd = NA,
          qaly = NA,
          qaly_sd = NA,
          icer = cycle_data$icer$avg
        )
      }
    }
  }
  
  # Convert the list to a data frame
  df <- bind_rows(data_list)
  return(df)
}

# 2. Prepare summary datasets
prepare_summary_data <- function(df) {
  # Time-series summary
  time_series <- df %>%
    filter(strategy != "icer") %>%
    group_by(cycle, strategy, state) %>%
    summarize(
      mean_count = mean(count, na.rm = TRUE),
      sd_count = sd(count, na.rm = TRUE),
      mean_total_cost = mean(total_cost, na.rm = TRUE),
      mean_surveillance_cost = mean(surveillance_cost, na.rm = TRUE),
      mean_qaly = mean(qaly, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Endpoint summary (cycle 20 or max cycle)
  max_cycle <- max(df$cycle)
  
  endpoint <- df %>%
    filter(strategy != "icer", cycle == max_cycle) %>%
    group_by(replication, strategy) %>%
    summarize(
      total_cost = first(total_cost),
      surveillance_cost = first(surveillance_cost),
      qaly = first(qaly),
      .groups = "drop"
    )
  
  # Cost summary by cycle
  cost_summary <- df %>%
    filter(strategy != "icer") %>%
    group_by(cycle, strategy) %>%
    summarize(
      total_cost = mean(total_cost, na.rm = TRUE),
      surveillance_cost = mean(surveillance_cost, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(total_cost, surveillance_cost),
      names_to = "cost_type",
      values_to = "cost"
    )
  
  return(list(
    time_series = time_series,
    endpoint = endpoint,
    cost_summary = cost_summary
  ))
}

# 3. Create plots
create_plots <- function(summary_data) {
  time_series <- summary_data$time_series
  endpoint <- summary_data$endpoint
  cost_summary <- summary_data$cost_summary
  
  # 1. Stacked area chart of mean state counts
  state_plot <- time_series %>%
    group_by(cycle, strategy) %>%
    mutate(total_count = sum(mean_count, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(proportion = mean_count / total_count) %>%
    ggplot(aes(x = cycle, y = mean_count, fill = state)) +
    geom_area(position = "stack") +
    facet_wrap(~ strategy) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "State Distribution over Time",
      x = "Cycle",
      y = "Number of Patients",
      fill = "Health State"
    ) +
    theme_minimal()
  
  # 2. Line plot of costs
  cost_plot <- ggplot(cost_summary, aes(x = cycle, y = cost, color = strategy, linetype = cost_type)) +
    geom_line(linewidth = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Mean Costs over Time",
      x = "Cycle",
      y = "Cost ($)",
      color = "Strategy",
      linetype = "Cost Type"
    ) +
    theme_minimal()
  
  # 3. Boxplots for endpoint costs and QALYs
  cost_boxplot <- ggplot(endpoint, aes(x = strategy, y = total_cost, fill = strategy)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Total Costs at Final Cycle",
      x = "Strategy",
      y = "Total Cost ($)",
      fill = "Strategy"
    ) +
    theme_minimal()
  
  qaly_boxplot <- ggplot(endpoint, aes(x = strategy, y = qaly, fill = strategy)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "QALYs at Final Cycle",
      x = "Strategy",
      y = "QALYs",
      fill = "Strategy"
    ) +
    theme_minimal()
  
  # 4. Scatter plot of cost vs. QALYs
  scatter_plot <- ggplot(endpoint, aes(x = qaly, y = total_cost, color = strategy)) +
    geom_point(alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Cost-Effectiveness Plane",
      x = "QALYs",
      y = "Total Cost ($)",
      color = "Strategy"
    ) +
    theme_minimal()
  
  # Return all plots
  return(list(
    state_plot = state_plot,
    cost_plot = cost_plot,
    cost_boxplot = cost_boxplot,
    qaly_boxplot = qaly_boxplot,
    scatter_plot = scatter_plot
  ))
}

# Main execution
main <- function() {
  # Find the most recent baseline result
  result_files <- list.files("results", pattern = "Baseline_Result_.*\\.json", full.names = TRUE)
  if (length(result_files) == 0) {
    stop("No baseline result files found in the results directory.")
  }
  
  # Sort files by modification time and get the most recent
  file_info <- file.info(result_files)
  most_recent <- result_files[which.max(file_info$mtime)]
  cat("Processing file:", most_recent, "\n")
  
  # Process the data
  df <- process_json_results(most_recent)
  summary_data <- prepare_summary_data(df)
  plots <- create_plots(summary_data)
  
  # Save the plots
  ggsave("results/state_distribution.png", plots$state_plot, width = 10, height = 6)
  ggsave("results/costs_over_time.png", plots$cost_plot, width = 10, height = 6)
  ggsave("results/cost_boxplot.png", plots$cost_boxplot, width = 8, height = 6)
  ggsave("results/qaly_boxplot.png", plots$qaly_boxplot, width = 8, height = 6)
  ggsave("results/cost_effectiveness_plane.png", plots$scatter_plot, width = 8, height = 6)
  
  # Save the processed data for further analysis
  saveRDS(df, "results/processed_data.rds")
  saveRDS(summary_data, "results/summary_data.rds")
  
  cat("Analysis complete. Plots and data saved to the results directory.\n")
  
  # Return invisible summary data
  invisible(summary_data)
}

# Run the main function if this script is executed directly
if (!interactive()) {
  main()
}