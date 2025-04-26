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
if (!require(future)) install.packages("future")
if (!require(promises)) install.packages("promises")
if (!require(shinyjs)) install.packages("shinyjs")

library(shiny)
library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(future)
library(promises)
library(shinyjs)

# Define health state names (needed for logging and visualizations)
v.n <- c("H1", "H2", "S1", "S2", "P", "D")

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
  useShinyjs(),
  titlePanel("Microsimulation Model GUI"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_reps", "Number of Replications:", 100, min = 1),
      numericInput("n_individuals", "Individuals per Simulation:", 1000, min = 1),
      numericInput("n_cycles", "Cycles per Simulation:", 20, min = 1),
      numericInput("time_limit", "Time Limit per Simulation (seconds):", 30, min = 1),
      numericInput("seed", "Random Seed:", 12345, min = 1),
      checkboxInput("verbose", "Verbose Output", FALSE),
      checkboxInput("sensitivity", "Run With Sensitivity Analysis", FALSE),
      actionButton("run_simulation", "Run Simulation", class = "btn-primary")
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

# Safe function to update reactive values from promises/futures
safeUpdate <- function(session, reactiveValue, value) {
  if (!is.null(session) && !is.null(reactiveValue)) {
    isolate({
      reactiveValue(value)
    })
    session$onFlushed(function() {
      # This ensures the UI is updated
    })
  }
}

server <- function(input, output, session) {
  # Initialize reactive values for storing simulation progress
  progress_logs <- reactiveVal("")
  results_summary <- reactiveVal("")
  current_plots_data <- reactiveVal(NULL)
  is_running <- reactiveVal(FALSE)

  # Safe function to append to progress log
  append_to_progress <- function(message) {
    # Get the current logs
    isolate({
      current <- progress_logs()
      # Update the reactive value
      progress_logs(paste0(current, message, "\n"))
    })
  }

  # Function to handle output from the model
  log_message <- function(message) {
    future::value(
      future({
        append_to_progress(message)
      })
    )
  }

  # Function to generate plots from json file
  generate_and_display_plots <- function(json_file) {
    tryCatch({
      # Process the data using functions from analyze_results.R
      df <- process_json_results(json_file)
      summary_data <- prepare_summary_data(df)
      
      # Store data for download handlers
      current_plots_data(summary_data)
      
      # Update the plots
      output$state_distribution <- renderPlot({
        req(summary_data$time_series)
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
        req(summary_data$cost_summary)
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
        req(summary_data$endpoint)
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
        req(summary_data$endpoint)
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
        req(summary_data$endpoint)
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
    }, error = function(e) {
      append_to_progress(paste("Error generating plots:", e$message))
    })
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

  # Custom print function that redirects to our log
  # This function will be called within the simulation process
  custom_cat <- function(...) {
    message <- paste0(...)
    session$sendCustomMessage(type = "log_message", message = message)
  }
  
  # Create a message handler to receive log messages
  observeEvent(input$log_message, {
    append_to_progress(input$log_message)
  })

  # Run simulation button handler
  observeEvent(input$run_simulation, {
    # Disable the button while running
    disable("run_simulation")
    is_running(TRUE)
    
    # Reset progress log and results summary
    progress_logs("")
    results_summary("")
    
    # Setup parameters
    n_reps <- input$n_reps
    n.i <- input$n_individuals
    n.t <- input$n_cycles
    verbose <- input$verbose
    seed <- input$seed
    sensitivity <- input$sensitivity
    time_limit <- if(sensitivity) input$time_limit else NULL
    
    # Create prefix for results file
    result_prefix <- if(sensitivity) "Sensitivity" else "Baseline"
    
    # Log simulation details
    append_to_progress(paste("Starting", result_prefix, "simulation..."))
    append_to_progress(paste("Replications:", n_reps))
    append_to_progress(paste("Individuals:", n.i))
    append_to_progress(paste("Cycles:", n.t))
    append_to_progress(paste("Sensitivity:", sensitivity))
    if (sensitivity) {
      append_to_progress(paste("Time limit per simulation:", time_limit, "seconds"))
    }
    append_to_progress(paste("Seed:", seed))
    append_to_progress("--------------------------------------")
    
    # Create a proxy output sink that will send messages to our log
    outputLog <- function(message) {
      isolate({
        current <- progress_logs()
        progress_logs(paste0(current, message, "\n"))
      })
      invisible(NULL)
    }
    
    # Create a function to update results incrementally
    updateIncrementalResults <- function(sim_index, sim_result) {
      if (!is.null(sim_result)) {
        msg <- sprintf("Simulation %d completed in %.2f seconds (ESD Cost: %.0f Surgery Cost: %.0f ESD QALYs: %.3f Surgery QALYs: %.3f ICER: %.0f )",
                      sim_index, 
                      sim_result$time_taken,
                      sim_result$esd$tc_hat,
                      sim_result$surgery$tc_hat,
                      sim_result$esd$te_hat,
                      sim_result$surgery$te_hat,
                      ifelse(is.finite(sim_result$icer), sim_result$icer, NA))
        
        # Add to progress log
        append_to_progress(msg)
        
        # Update the incremental results
        isolate({
          current <- results_summary()
          # If this is first simulation, start with a header
          if (current == "") {
            current <- paste0(
              result_prefix, " Simulation Incremental Results:\n",
              "--------------------------------\n"
            )
          }
          results_summary(paste0(current, "Sim ", sim_index, ": ESD Cost: ", round(sim_result$esd$tc_hat), 
                               ", Surgery Cost: ", round(sim_result$surgery$tc_hat),
                               ", ICER: ", round(sim_result$icer), "\n"))
        })
      }
    }
    
    # Create a future to run the simulation
    simulation_future <- future({
      # Define v.n directly inside the future instead of trying to capture it from parent environment
      v.n <- c("H1", "H2", "S1", "S2", "P", "D")
      
      # Capture cat output
      con <- textConnection("output_capture", "w", local = TRUE)
      sink(con, append = TRUE)
      
      # Run the simulation
      result <- run_multiple_simulations(
        n_reps = n_reps,
        sensitivity = sensitivity,
        time_limit_seconds = time_limit,
        n.i = n.i, 
        n.t = n.t, 
        verbose = verbose, 
        seed = seed
      )
      
      # Stop capturing output
      sink()
      close(con)
      
      # Return both the simulation result and captured output
      list(result = result, output = output_capture, v.n = v.n)
    }, seed = TRUE) # Add seed=TRUE to properly handle random numbers in the future
    
    # Handle the future when it completes
    simulation_future %...>% 
      (function(result_data) {
        # Process the captured output
        for (line in result_data$output) {
          isolate({
            current <- progress_logs()
            progress_logs(paste0(current, line, "\n"))
          })
        }
        
        # Extract the simulation result and captured v.n
        result <- result_data$result
        v.n <- result_data$v.n
        
        # Log completion
        isolate({
          current <- progress_logs()
          progress_logs(paste0(current, "--------------------------------------\n",
                              "Completed ", result$completed_sims, " of ", n_reps, " simulations!\n"))
        })
        
        # Prepare json file
        json_file <- get_timestamped_filename(result_prefix)
        
        # Check if we have completed simulations
        if (result$completed_sims > 0) {
          # Create compatible structure for the logging function
          log_results <- list(all_results = list())
          
          # Get completed simulations only
          for (i in 1:min(n_reps, result$completed_sims)) {
            if (!is.null(result$results[[i]]) && result$results[[i]]$completed) {
              log_results$all_results[[i]] <- result$results[[i]]
            }
          }
          
          # Only log if we have valid results
          if (length(log_results$all_results) > 0) {
            log_simulation_cycle_results(log_results, n.i, v.n, json_file)
            isolate({
              current <- progress_logs()
              progress_logs(paste0(current, "Results saved to: ", json_file, "\n"))
            })
            
            # Generate and display plots
            generate_and_display_plots(json_file)
          }
        }
        
        # Set final summary results
        if (result$completed_sims > 0 && !is.null(result$summary)) {
          summary <- result$summary
          final_summary <- paste0(
            result_prefix, " Simulation Results:\n",
            "Completed Simulations: ", result$completed_sims, " of ", n_reps, "\n",
            "Using seed: ", seed, "\n",
            "ESD Mean QALYs: ", round(summary$mean_esd_qaly, 4), 
              " (SD: ", round(summary$sd_esd_qaly, 4), ")\n",
            "Surgery Mean QALYs: ", round(summary$mean_surgery_qaly, 4), 
              " (SD: ", round(summary$sd_surgery_qaly, 4), ")\n",
            "ESD Mean Cost: ", round(summary$mean_esd_cost, 2), 
              " (SD: ", round(summary$sd_esd_cost, 2), ")\n",
            "Surgery Mean Cost: ", round(summary$mean_surgery_cost, 2), 
              " (SD: ", round(summary$sd_surgery_cost, 2), ")\n",
            "Delta QALYs: ", round(summary$mean_delta_qaly, 4), 
              " (SD: ", round(summary$sd_delta_qaly, 4), ")\n",
            "Delta Cost: ", round(summary$mean_delta_cost, 2), 
              " (SD: ", round(summary$sd_delta_cost, 2), ")\n",
            "ICER: ", round(summary$mean_icer, 2), 
              " (SD: ", round(summary$sd_icer, 2), ")\n"
          )
          
          # Add percentiles if available
          if (!is.null(summary$icer_percentiles)) {
            final_summary <- paste0(
              final_summary,
              "\nICER Percentiles:\n",
              "  2.5%: ", round(summary$icer_percentiles[1], 2), "\n",
              "  25%: ", round(summary$icer_percentiles[2], 2), "\n",
              "  50% (median): ", round(summary$icer_percentiles[3], 2), "\n",
              "  75%: ", round(summary$icer_percentiles[4], 2), "\n",
              "  97.5%: ", round(summary$icer_percentiles[5], 2), "\n"
            )
          }
          
          if (length(log_results$all_results) > 0) {
            final_summary <- paste0(final_summary, "Cycle log saved to: ", json_file, "\n")
          }
          
          results_summary(final_summary)
        } else {
          results_summary(paste0(
            "No simulations completed successfully. ",
            if(sensitivity) "Try increasing the time limit." else "", 
            "\n"
          ))
        }
      }) %...!% 
      (function(error) {
        # Handle any errors in the future
        isolate({
          current <- progress_logs()
          progress_logs(paste0(current, "Error: ", as.character(error), "\n"))
        })
      }) %...>%
      (function(...) {
        # Always re-enable the button when done
        enable("run_simulation")
        is_running(FALSE)
      })
  })
  
  # Configure future to run in a separate R session
  plan(multisession)
  
  # Setup custom message handler
  session$registerDataObj(
    name = "log_message",
    data = {
      function(message, id = NULL) {
        isolate({
          current <- progress_logs()
          progress_logs(paste0(current, message, "\n"))
        })
      }
    },
    filterFunc = function(message, id = NULL) {
      return(TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
