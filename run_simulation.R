# run_simulation.R
# Script to execute simulation runs and post-processing (separated from model definition)

# Source the model and required functions
source("./microsimulation_model.R")
source("./functions.R")

# Shiny GUI for running baseline and sensitivity analysis
if (!require(shiny)) install.packages("shiny"); library(shiny)
if (!require(jsonlite)) install.packages("jsonlite"); library(jsonlite)

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
      actionButton("run_baseline", "Run Baseline Simulation"),
      actionButton("run_sensitivity", "Run Sensitivity Analysis")
    ),
    mainPanel(
      verbatimTextOutput("sim_output")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$run_baseline, {
    n_reps <- input$n_reps
    n.i <- input$n_individuals
    n.t <- input$n_cycles
    verbose <- input$verbose
    seed <- input$seed
    
    # Pass the seed to run_baseline_simulation
    result <- run_baseline_simulation(
      n.i = n.i, 
      n.t = n.t, 
      verbose = verbose, 
      n_reps = n_reps,
      seed = seed
    )
    
    # Logging for simulation results
    json_file <- get_timestamped_filename("Baseline")
    log_simulation_cycle_results(result, n.i, v.n, json_file)
    
    if (n_reps == 1) {
      output$sim_output <- renderPrint({
        cat("Baseline Simulation Results (Single Run):\n")
        cat("Using seed:", seed, "\n")
        cat("ESD Mean QALYs:", mean(result$esd$te), "\n")
        cat("Surgery Mean QALYs:", mean(result$surgery$te), "\n")
        cat("ESD Mean Cost:", mean(result$esd$tc), "\n")
        cat("Surgery Mean Cost:", mean(result$surgery$tc), "\n")
        cat("Delta QALYs:", result$delta_e, "\n")
        cat("Delta Cost:", result$delta_c, "\n")
        cat("ICER:", result$icer, "\n")
        cat("Cycle log saved to:", json_file, "\n")
      })
    } else {
      output$sim_output <- renderPrint({
        cat("Baseline Simulation Replications:", result$n_reps, "\n")
        cat("Using seed:", seed, "\n")
        cat("ESD Mean QALYs:", result$esd_qalys_mean, "(SD:", result$esd_qalys_sd, ")\n")
        cat("Surgery Mean QALYs:", result$surgery_qalys_mean, "(SD:", result$surgery_qalys_sd, ")\n")
        cat("ESD Mean Cost:", result$esd_costs_mean, "(SD:", result$esd_costs_sd, ")\n")
        cat("Surgery Mean Cost:", result$surgery_costs_mean, "(SD:", result$surgery_costs_sd, ")\n")
        cat("Delta QALYs:", result$delta_qalys_mean, "(SD:", result$delta_qalys_sd, ")\n")
        cat("Delta Cost:", result$delta_costs_mean, "(SD:", result$delta_costs_sd, ")\n")
        cat("ICER:", result$icer_mean, "(SD:", result$icer_sd, ")\n")
        cat("Cycle log saved to:", json_file, "\n")
      })
    }
  })

  observeEvent(input$run_sensitivity, {
    seed <- input$seed
    n_reps <- input$n_reps
    
    result <- sensitivity_analysis(
      n_reps = n_reps,
      time_limit_seconds = input$time_limit,
      n.i = input$n_individuals,
      n.t = input$n_cycles,
      verbose = input$verbose,
      seed = seed
    )
    
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
      }
    }
    
    output$sim_output <- renderPrint({
      cat("Sensitivity Analysis Results:\n")
      cat("Using seed:", seed, "\n")
      cat("Completed Simulations:", result$completed_sims, "\n")
      
      if (result$completed_sims > 0) {
        cat("Mean ICER:", result$mean_icer, "\n")
        cat("Mean Incremental Cost:", result$mean_delta_cost, "\n")
        cat("Mean Incremental QALYs:", result$mean_delta_qaly, "\n")
        cat("ICER Percentiles:\n")
        print(result$icer_percentiles)
        if (length(sensitivity_results$all_results) > 0) {
          cat("Cycle log saved to:", json_file, "\n") 
        } else {
          cat("No complete results available for cycle logging\n")
        }
      } else {
        cat("No simulations completed successfully. Try increasing the time limit.\n")
      }
    })
  })
}

shinyApp(ui = ui, server = server)