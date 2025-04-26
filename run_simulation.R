# run_simulation.R
# Script to execute simulation runs and post-processing (separated from model definition)

# Source the model and required functions
source("./microsimulation_model.R")
source("./functions.R")

# Install and load required packages
if (!require(shiny)) install.packages("shiny")
if (!require(jsonlite)) install.packages("jsonlite")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(ggplot2)) install.packages("ggplot2")

library(shiny)
library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)

# Process JSON results function (moved from analyze_results.R)
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

# Prepare summary datasets function (moved from analyze_results.R)
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

get_timestamped_filename <- function(prefix) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_dir <- "./results"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  file.path(results_dir, paste0(prefix, "_Result_", ts, ".json"))
}

ui <- fluidPage(
  titlePanel("Microsimulation Model GUI"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_reps", "Number of Replications:", 100, min = 1),
      numericInput("n_individuals", "Individuals per Simulation:", 1000, min = 1),
      numericInput("n_cycles", "Cycles per Simulation:", 20, min = 1),
      numericInput("time_limit", "Time Limit per Simulation (seconds):", 30, min = 1),
      numericInput("seed", "Random Seed:", 12345, min = 1),
      checkboxInput("verbose", "Verbose Output", FALSE),
      actionButton("run_baseline", "Run Baseline Simulation", class = "btn-primary"),
      actionButton("run_sensitivity", "Run Sensitivity Analysis", class = "btn-info")
    ),
    mainPanel(
      tabsetPanel(
        id = "result_tabs",
        tabPanel(
          "Results Summary", 
          # Fixed size scrollable output for simulation progress
          tags$div(
            style = "height:300px; overflow-y: scroll; background-color: #f5f5f5; padding: 10px; border: 1px solid #ddd; margin-bottom: 20px;",
            verbatimTextOutput("sim_progress")
          ),
          # Final results summary in a separate panel
          tags$div(
            style = "background-color: #fff; padding: 15px; border: 1px solid #ddd;",
            h4("Simulation Results Summary"),
            verbatimTextOutput("sim_results")
          )
        ),
        tabPanel(
          "Visualizations",
          h3("Simulation Visualizations"),
          
          # State distribution plot with download button
          fluidRow(
            column(12, 
              tags$div(class = "plot-container", style = "position: relative;",
                plotOutput("state_distribution", height = "400px"),
                tags$div(class = "download-btn", style = "position: absolute; top: 10px; right: 10px; z-index: 100;",
                  downloadButton("download_state_dist", "", icon = icon("download"), 
                                style = "padding: 5px; background-color: rgba(255, 255, 255, 0.7);")
                )
              )
            )
          ),
          
          # Costs plot with download button
          fluidRow(
            column(12,
              tags$div(class = "plot-container", style = "position: relative;",
                plotOutput("costs_plot", height = "400px"),
                tags$div(class = "download-btn", style = "position: absolute; top: 10px; right: 10px; z-index: 100;",
                  downloadButton("download_costs", "", icon = icon("download"),
                                style = "padding: 5px; background-color: rgba(255, 255, 255, 0.7);")
                )
              )
            )
          ),
          
          # Cost and QALY boxplots with download buttons
          fluidRow(
            column(6, 
              tags$div(class = "plot-container", style = "position: relative;",
                plotOutput("cost_boxplot", height = "350px"),
                tags$div(class = "download-btn", style = "position: absolute; top: 10px; right: 10px; z-index: 100;",
                  downloadButton("download_cost_box", "", icon = icon("download"),
                               style = "padding: 5px; background-color: rgba(255, 255, 255, 0.7);")
                )
              )
            ),
            column(6, 
              tags$div(class = "plot-container", style = "position: relative;",
                plotOutput("qaly_boxplot", height = "350px"),
                tags$div(class = "download-btn", style = "position: absolute; top: 10px; right: 10px; z-index: 100;",
                  downloadButton("download_qaly_box", "", icon = icon("download"),
                               style = "padding: 5px; background-color: rgba(255, 255, 255, 0.7);")
                )
              )
            )
          ),
          
          # Scatter plot with download button
          fluidRow(
            column(12,
              tags$div(class = "plot-container", style = "position: relative;",
                plotOutput("scatter_plot", height = "400px"),
                tags$div(class = "download-btn", style = "position: absolute; top: 10px; right: 10px; z-index: 100;",
                  downloadButton("download_scatter", "", icon = icon("download"),
                               style = "padding: 5px; background-color: rgba(255, 255, 255, 0.7);")
                )
              )
            )
          ),
          
          # CSS to show download button only on hover
          tags$style(HTML("
            .plot-container .download-btn { 
              opacity: 0;
              transition: opacity 0.3s;
            }
            .plot-container:hover .download-btn { 
              opacity: 1; 
            }
          "))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Initialize reactive values for storing simulation progress
  progress_logs <- reactiveVal("")
  results_summary <- reactiveVal("")
  current_plots_data <- reactiveVal(NULL)

  # Function to append to progress log
  append_to_progress <- function(message) {
    current <- progress_logs()
    progress_logs(paste0(current, message, "\n"))
  }

  # Function to generate plots from json file
  generate_and_display_plots <- function(json_file) {
    # Process the data using functions from analyze_results.R
    df <- process_json_results(json_file)
    summary_data <- prepare_summary_data(df)
    
    # Store data for download handlers
    current_plots_data(summary_data)
    
    # Update the plots
    output$state_distribution <- renderPlot({
      summary_data$time_series %>%
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
    })
    
    output$costs_plot <- renderPlot({
      summary_data$cost_summary %>%
        ggplot(aes(x = cycle, y = cost, color = strategy, linetype = cost_type)) +
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
    })
    
    output$cost_boxplot <- renderPlot({
      summary_data$endpoint %>%
        ggplot(aes(x = strategy, y = total_cost, fill = strategy)) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Set1") +
        labs(
          title = "Total Costs at Final Cycle",
          x = "Strategy",
          y = "Total Cost ($)",
          fill = "Strategy"
        ) +
        theme_minimal()
    })
    
    output$qaly_boxplot <- renderPlot({
      summary_data$endpoint %>%
        ggplot(aes(x = strategy, y = qaly, fill = strategy)) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Set1") +
        labs(
          title = "QALYs at Final Cycle",
          x = "Strategy",
          y = "QALYs",
          fill = "Strategy"
        ) +
        theme_minimal()
    })
    
    output$scatter_plot <- renderPlot({
      summary_data$endpoint %>%
        ggplot(aes(x = qaly, y = total_cost, color = strategy)) +
        geom_point(alpha = 0.7) +
        scale_color_brewer(palette = "Set1") +
        labs(
          title = "Cost-Effectiveness Plane",
          x = "QALYs",
          y = "Total Cost ($)",
          color = "Strategy"
        ) +
        theme_minimal()
    })
    
    # Switch to the visualizations tab
    updateTabsetPanel(session, "result_tabs", selected = "Visualizations")
  }
  
  # Connect the output objects to the reactive values
  output$sim_progress <- renderText({
    progress_logs()
  })
  
  output$sim_results <- renderText({
    results_summary()
  })
  
  # Download handlers for each plot
  output$download_state_dist <- downloadHandler(
    filename = function() {
      paste("state_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      summary_data <- current_plots_data()
      if (!is.null(summary_data)) {
        ggsave(file, plot = {
          summary_data$time_series %>%
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
        }, width = 10, height = 7, dpi = 300)
      }
    }
  )
  
  output$download_costs <- downloadHandler(
    filename = function() {
      paste("costs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      summary_data <- current_plots_data()
      if (!is.null(summary_data)) {
        ggsave(file, plot = {
          summary_data$cost_summary %>%
            ggplot(aes(x = cycle, y = cost, color = strategy, linetype = cost_type)) +
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
        }, width = 10, height = 7, dpi = 300)
      }
    }
  )
  
  output$download_cost_box <- downloadHandler(
    filename = function() {
      paste("cost_boxplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      summary_data <- current_plots_data()
      if (!is.null(summary_data)) {
        ggsave(file, plot = {
          summary_data$endpoint %>%
            ggplot(aes(x = strategy, y = total_cost, fill = strategy)) +
            geom_boxplot() +
            scale_fill_brewer(palette = "Set1") +
            labs(
              title = "Total Costs at Final Cycle",
              x = "Strategy",
              y = "Total Cost ($)",
              fill = "Strategy"
            ) +
            theme_minimal()
        }, width = 8, height = 6, dpi = 300)
      }
    }
  )
  
  output$download_qaly_box <- downloadHandler(
    filename = function() {
      paste("qaly_boxplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      summary_data <- current_plots_data()
      if (!is.null(summary_data)) {
        ggsave(file, plot = {
          summary_data$endpoint %>%
            ggplot(aes(x = strategy, y = qaly, fill = strategy)) +
            geom_boxplot() +
            scale_fill_brewer(palette = "Set1") +
            labs(
              title = "QALYs at Final Cycle",
              x = "Strategy",
              y = "QALYs",
              fill = "Strategy"
            ) +
            theme_minimal()
        }, width = 8, height = 6, dpi = 300)
      }
    }
  )
  
  output$download_scatter <- downloadHandler(
    filename = function() {
      paste("cost_effectiveness_plane_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png", sep = "")
    },
    content = function(file) {
      summary_data <- current_plots_data()
      if (!is.null(summary_data)) {
        ggsave(file, plot = {
          summary_data$endpoint %>%
            ggplot(aes(x = qaly, y = total_cost, color = strategy)) +
            geom_point(alpha = 0.7) +
            scale_color_brewer(palette = "Set1") +
            labs(
              title = "Cost-Effectiveness Plane",
              x = "QALYs",
              y = "Total Cost ($)",
              color = "Strategy"
            ) +
            theme_minimal()
        }, width = 10, height = 7, dpi = 300)
      }
    }
  )

  # Simulate baseline model
  observeEvent(input$run_baseline, {
    # Reset progress log
    progress_logs("")
    append_to_progress("Starting baseline simulation...")
    append_to_progress(paste("Replications:", input$n_reps))
    append_to_progress(paste("Individuals:", input$n_individuals))
    append_to_progress(paste("Cycles:", input$n_cycles))
    append_to_progress(paste("Seed:", input$seed))
    append_to_progress("--------------------------------------")
    
    # Setup parameters
    n_reps <- input$n_reps
    n.i <- input$n_individuals
    n.t <- input$n_cycles
    verbose <- input$verbose
    seed <- input$seed
    
    # Create progress tracker that updates UI
    progress_tracker <- function(i, completed = FALSE) {
      if (completed) {
        append_to_progress(paste("Simulation", i, "of", n_reps, "completed."))
      } else {
        append_to_progress(paste("Running simulation", i, "of", n_reps, "..."))
      }
    }
    
    # Run the simulation with progress tracking
    withProgress(message = "Running baseline simulation", {
      for (i in 1:n_reps) {
        incProgress(1/n_reps, detail = paste("Simulation", i, "of", n_reps))
        progress_tracker(i)
      }
      
      # Actually run the simulation
      result <- run_baseline_simulation(
        n.i = n.i, 
        n.t = n.t, 
        verbose = verbose, 
        n_reps = n_reps,
        seed = seed
      )
      
      # Log completion
      append_to_progress("--------------------------------------")
      append_to_progress("Simulation completed!")
      
      # Logging for simulation results
      json_file <- get_timestamped_filename("Baseline")
      log_simulation_cycle_results(result, n.i, v.n, json_file)
      append_to_progress(paste("Results saved to:", json_file))
      
      # Set final summary results
      if (n_reps == 1) {
        final_summary <- paste0(
          "Baseline Simulation Results (Single Run):\n",
          "Using seed: ", seed, "\n",
          "ESD Mean QALYs: ", round(mean(result$esd$te), 4), "\n",
          "Surgery Mean QALYs: ", round(mean(result$surgery$te), 4), "\n",
          "ESD Mean Cost: ", round(mean(result$esd$tc), 2), "\n",
          "Surgery Mean Cost: ", round(mean(result$surgery$tc), 2), "\n",
          "Delta QALYs: ", round(result$delta_e, 4), "\n",
          "Delta Cost: ", round(result$delta_c, 2), "\n",
          "ICER: ", round(result$icer, 2), "\n",
          "Cycle log saved to: ", json_file, "\n"
        )
      } else {
        final_summary <- paste0(
          "Baseline Simulation Replications: ", result$n_reps, "\n",
          "Using seed: ", seed, "\n",
          "ESD Mean QALYs: ", round(result$esd_qalys_mean, 4), " (SD: ", round(result$esd_qalys_sd, 4), ")\n",
          "Surgery Mean QALYs: ", round(result$surgery_qalys_mean, 4), " (SD: ", round(result$surgery_qalys_sd, 4), ")\n",
          "ESD Mean Cost: ", round(result$esd_costs_mean, 2), " (SD: ", round(result$esd_costs_sd, 2), ")\n",
          "Surgery Mean Cost: ", round(result$surgery_costs_mean, 2), " (SD: ", round(result$surgery_costs_sd, 2), ")\n",
          "Delta QALYs: ", round(result$delta_qalys_mean, 4), " (SD: ", round(result$delta_qalys_sd, 4), ")\n",
          "Delta Cost: ", round(result$delta_costs_mean, 2), " (SD: ", round(result$delta_costs_sd, 2), ")\n",
          "ICER: ", round(result$icer_mean, 2), " (SD: ", round(result$icer_sd, 2), ")\n",
          "Cycle log saved to: ", json_file, "\n"
        )
      }
      results_summary(final_summary)
      
      # Generate and display plots for the baseline simulation
      generate_and_display_plots(json_file)
    })
  })

  # Run sensitivity analysis
  observeEvent(input$run_sensitivity, {
    # Reset progress log
    progress_logs("")
    append_to_progress("Starting sensitivity analysis...")
    append_to_progress(paste("Replications:", input$n_reps))
    append_to_progress(paste("Individuals:", input$n_individuals))
    append_to_progress(paste("Cycles:", input$n_cycles))
    append_to_progress(paste("Time limit per simulation:", input$time_limit, "seconds"))
    append_to_progress(paste("Seed:", input$seed))
    append_to_progress("--------------------------------------")
    
    # Setup parameters
    seed <- input$seed
    n_reps <- input$n_reps
    time_limit <- input$time_limit
    
    # Create a sensitivity analysis with progress tracking
    withProgress(message = "Running sensitivity analysis", {
      # Run the sensitivity analysis with progress updates
      result <- sensitivity_analysis(
        n_reps = n_reps,
        time_limit_seconds = time_limit,
        n.i = input$n_individuals,
        n.t = input$n_cycles,
        verbose = input$verbose,
        seed = seed
      )
      
      # Every few seconds, check progress
      for (i in 1:n_reps) {
        incProgress(1/n_reps, detail = paste("Simulation", i, "of", n_reps))
        if (i %% 5 == 0 || i == n_reps) {
          append_to_progress(paste("Completed", i, "of", n_reps, "simulations..."))
        }
        Sys.sleep(0.1)  # Small delay to allow UI updates
      }
      
      # Log completion
      append_to_progress("--------------------------------------")
      append_to_progress(paste("Sensitivity analysis completed with", result$completed_sims, "of", n_reps, "simulations!"))
      
      # Logging for sensitivity analysis results - handle different structure
      json_file <- get_timestamped_filename("Sensitivity")
      
      # Check if we have completed simulations to process
      if (result$completed_sims > 0) {
        # For sensitivity analysis, results are stored differently
        # Create compatible structure for the logging function
        sensitivity_results <- list(all_results = list())
        
        # Get completed simulations only
        for (i in 1:min(n_reps, result$completed_sims)) {
          if (!is.null(result$results[[i]]) && result$results[[i]]$completed) {
            sensitivity_results$all_results[[i]] <- result$results[[i]]
          }
        }
        
        # Only log if we have valid results
        if (length(sensitivity_results$all_results) > 0) {
          log_simulation_cycle_results(sensitivity_results, input$n_individuals, v.n, json_file)
          append_to_progress(paste("Results saved to:", json_file))
          
          # Generate and display plots for the sensitivity analysis
          generate_and_display_plots(json_file)
        }
      }
      
      # Set final summary results
      final_summary <- paste0(
        "Sensitivity Analysis Results:\n",
        "Using seed: ", seed, "\n",
        "Completed Simulations: ", result$completed_sims, "\n"
      )
      
      if (result$completed_sims > 0) {
        final_summary <- paste0(
          final_summary,
          "Mean ICER: ", round(result$mean_icer, 2), "\n",
          "Mean Incremental Cost: ", round(result$mean_delta_cost, 2), "\n",
          "Mean Incremental QALYs: ", round(result$mean_delta_qaly, 4), "\n",
          "ICER Percentiles:\n",
          "  2.5%: ", round(result$icer_percentiles[1], 2), "\n",
          "  25%: ", round(result$icer_percentiles[2], 2), "\n",
          "  50% (median): ", round(result$icer_percentiles[3], 2), "\n",
          "  75%: ", round(result$icer_percentiles[4], 2), "\n",
          "  97.5%: ", round(result$icer_percentiles[5], 2), "\n"
        )
        
        if (length(sensitivity_results$all_results) > 0) {
          final_summary <- paste0(final_summary, "Cycle log saved to: ", json_file, "\n")
        } else {
          final_summary <- paste0(final_summary, "No complete results available for cycle logging\n")
        }
      } else {
        final_summary <- paste0(final_summary, "No simulations completed successfully. Try increasing the time limit.\n")
      }
      
      results_summary(final_summary)
    })
  })
}

shinyApp(ui = ui, server = server)
