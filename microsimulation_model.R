############################################################################################
################# Microsimulation modeling using R: a tutorial #### 2018 ###################
############################################################################################
# This code forms the basis for the microsimulation model of the article: 
#
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22.
#
# Please cite the article when using this code
# 
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial
#
# To program this tutorial we made use of 
# R: 3.3.0 GUI 1.68 Mavericks build (7202)
# RStudio: Version 1.0.136 2009-2016 RStudio, Inc.

############################################################################################
################# Code of Appendix A #######################################################
############################################################################################
# Clear environment variables
rm(list = ls())  # remove any variables in R's memory 

# First, source all required files BEFORE defining functions
# Source the files with exact case matching to avoid path issues
source("./functions.r")  # lowercase 'r' to match actual filename
source("./BetaParmsFromQuantiles.R")
source("./GammaParmsFromQuantiles.R")

##################################### Function to run a timed simulation ###################
run_timed_simulation <- function(time_limit_seconds = 60, n.i = 1000, n.t = 30, verbose = FALSE, sim_index = 1) {
  # Start timing
  start_time <- Sys.time()
  
  # Initialize model parameters 
  initialize_model_parameters(force = TRUE)
  
  # Use simulation index to create a unique seed for each simulation
  set.seed(12345 + sim_index)
  
  # Generate random values for probabilities, costs, and utilities
  generate_new_random_probs(print_output = FALSE)
  generate_new_random_c_values(print_output = FALSE)
  generate_new_random_u_values(print_output = FALSE)

  
  # Define the initial state vectors
  v.M_1 <- rep("H1", n.i)  # all start in the H1 (ESD) state
  v.M_2 <- rep("H2", n.i)  # all start in the H2 (Surgery) state
  
  # Treatment names
  v.Trt <- c("ESD", "Surgery")
  
  # Discount rates
  d.c <- d.e <- 0.03  # discount rate for costs and QALYs
  
  # Run ESD and Surgery simulations
  if (verbose) cat("Running ESD simulation:\n")
  sim_esd <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = (12345 + sim_index*2))
  
  # Check if we've exceeded the time limit
  if (difftime(Sys.time(), start_time, units = "secs") > time_limit_seconds) {
    if (verbose) cat("Time limit reached, returning partial results\n")
    return(list(
      esd = sim_esd,
      surgery = NULL,
      completed = FALSE,
      runtime = difftime(Sys.time(), start_time, units = "secs")
    ))
  }
  
  if (verbose) cat("Running Surgery simulation:\n")
  sim_surgery <- MicroSim(v.M_2, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = (12345 + sim_index*2 + 1))
  
  # Calculate incremental values
  delta_c <- mean(sim_surgery$tc) - mean(sim_esd$tc)
  delta_e <- mean(sim_surgery$te) - mean(sim_esd$te)
  icer <- delta_c / delta_e
  
  # Return the simulation results
  return(list(
    esd = sim_esd,
    surgery = sim_surgery,
    delta_c = delta_c,
    delta_e = delta_e,
    icer = icer,
    completed = TRUE,
    runtime = difftime(Sys.time(), start_time, units = "secs")
  ))
}

##################################### Run multiple simulations ##############################
run_multiple_simulations <- function(n_sims = 1000, time_limit_seconds = 60, n.i = 1000, n.t = 30, verbose = FALSE) {
  # Prepare storage for results
  results <- list()
  completed_sims <- 0
  total_start_time <- Sys.time()
  
  cat("Starting", n_sims, "simulations with time limit of", time_limit_seconds, "seconds each\n")
  cat("Total individuals per simulation:", n.i, "\n")
  cat("Total cycles per simulation:", n.t, "\n\n")
  
  # Create results directory if it doesn't exist
  results_dir <- "./results"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  
  # Initialize results file with headers
  results_file <- file.path(results_dir, "simulation_results_live.csv")
  write("simulation,esd_qalys,surgery_qalys,esd_cost,surgery_cost,delta_qalys,delta_cost,icer", 
        file = results_file)
  
  for (i in 1:n_sims) {
    cat("Starting simulation", i, "of", n_sims, "\n")
    sim_start_time <- Sys.time()
    
    # Run a single timed simulation
    sim_result <- run_timed_simulation(time_limit_seconds, n.i, n.t, verbose, sim_index = i)
    
    # Store results
    results[[i]] <- sim_result
    
    # Update completion counter
    if (sim_result$completed) {
      completed_sims <- completed_sims + 1
    }
    
    # Print progress with more detailed information
    sim_runtime <- difftime(Sys.time(), sim_start_time, units = "secs")
    cat("Simulation", i, "completed in", round(sim_runtime, 2), "seconds")
    
    if (sim_result$completed) {
      # Calculate mean costs and QALYs for reporting
      esd_qalys <- mean(sim_result$esd$te)
      surgery_qalys <- mean(sim_result$surgery$te)
      esd_cost <- mean(sim_result$esd$tc)
      surgery_cost <- mean(sim_result$surgery$tc)
      
      # Print detailed output including costs
      cat(" (ESD Cost:", round(esd_cost, 0), 
          "Surgery Cost:", round(surgery_cost, 0),
          "ESD QALYs:", round(esd_qalys, 3), 
          "Surgery QALYs:", round(surgery_qalys, 3), 
          "ICER:", round(sim_result$icer, 0), ")\n")
      
      # Write results to file line by line immediately after each simulation
      result_line <- paste(i, esd_qalys, surgery_qalys, esd_cost, surgery_cost, 
                          sim_result$delta_e, sim_result$delta_c, sim_result$icer, 
                          sep = ",")
      write(result_line, file = results_file, append = TRUE)
    } else {
      cat(" (incomplete - time limit reached)\n")
    }
  }
  
  total_runtime <- difftime(Sys.time(), total_start_time, units = "secs")
  cat("\nAll simulations completed in", round(total_runtime, 2), "seconds\n")
  cat("Completed simulations:", completed_sims, "out of", n_sims, "\n\n")
  
  # Compile summary statistics
  icers <- sapply(results[1:completed_sims], function(x) if(x$completed) x$icer else NA)
  delta_costs <- sapply(results[1:completed_sims], function(x) if(x$completed) x$delta_c else NA)
  delta_qalys <- sapply(results[1:completed_sims], function(x) if(x$completed) x$delta_e else NA)
  
  # Remove NA values
  icers <- icers[!is.na(icers)]
  delta_costs <- delta_costs[!is.na(delta_costs)]
  delta_qalys <- delta_qalys[!is.na(delta_qalys)]
  
  # Print summary statistics
  cat("Summary of completed simulations:\n")
  cat("Mean ICER:", round(mean(icers), 2), "(SD:", round(sd(icers), 2), ")\n")
  cat("Mean Incremental Cost:", round(mean(delta_costs), 2), "(SD:", round(sd(delta_costs), 2), ")\n")
  cat("Mean Incremental QALYs:", round(mean(delta_qalys), 4), "(SD:", round(sd(delta_qalys), 4), ")\n")
  
  # Percentiles for ICER
  icer_percentiles <- quantile(icers, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  cat("\nICER Percentiles:\n")
  cat("2.5%:", round(icer_percentiles[1], 2), "\n")
  cat("25%:", round(icer_percentiles[2], 2), "\n")
  cat("50% (median):", round(icer_percentiles[3], 2), "\n")
  cat("75%:", round(icer_percentiles[4], 2), "\n")
  cat("97.5%:", round(icer_percentiles[5], 2), "\n")
  
  return(list(
    results = results,
    completed_sims = completed_sims,
    total_runtime = total_runtime,
    mean_icer = mean(icers),
    sd_icer = sd(icers),
    mean_delta_cost = mean(delta_costs),
    sd_delta_cost = sd(delta_costs),
    mean_delta_qaly = mean(delta_qalys),
    sd_delta_qaly = sd(delta_qalys),
    icer_percentiles = icer_percentiles
  ))
}

# Define functions first, THEN run the simulation
# ===============================================================================

# Make all objects available in global environment before continuing
assign("run_timed_simulation", run_timed_simulation, envir = .GlobalEnv)
assign("run_multiple_simulations", run_multiple_simulations, envir = .GlobalEnv)

# Now execute the simulation
# ===============================================================================

# Set simulation parameters
n_simulations <- 100          # Number of simulations to run
time_limit_per_sim <- 30      # Time limit per simulation in seconds
individuals_per_sim <- 1000   # Number of individuals per simulation
cycles_per_sim <- 30          # Number of cycles per simulation
verbose_output <- FALSE       # Whether to show verbose output

cat("Starting simulation with the following parameters:\n")
cat("Number of simulations:", n_simulations, "\n")
cat("Time limit per simulation:", time_limit_per_sim, "seconds\n")
cat("Individuals per simulation:", individuals_per_sim, "\n")
cat("Cycles per simulation:", cycles_per_sim, "\n")

# Run the multiple simulations
multi_sim_results <- run_multiple_simulations(
  n_sims = n_simulations,
  time_limit_seconds = time_limit_per_sim,
  n.i = individuals_per_sim,
  n.t = cycles_per_sim,
  verbose = verbose_output
)

# Optional: Save results to file
results_dir <- "./results"
# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
save(multi_sim_results, file = file.path(results_dir, "simulation_results.RData"))

# Enhanced post-processing to store detailed results
# Extract and save all QALYs and ICERs
all_results <- data.frame(
  simulation = integer(),
  esd_qalys_mean = numeric(),
  esd_qalys_sd = numeric(),
  surgery_qalys_mean = numeric(), 
  surgery_qalys_sd = numeric(),
  esd_costs_mean = numeric(),
  esd_costs_sd = numeric(),
  surgery_costs_mean = numeric(),
  surgery_costs_sd = numeric(),
  delta_qalys = numeric(),
  delta_costs = numeric(),
  icer = numeric()
)

# Populate the data frame with results from each simulation
for (i in 1:multi_sim_results$completed_sims) {
  if (multi_sim_results$results[[i]]$completed) {
    all_results[i, "simulation"] <- i
    all_results[i, "esd_qalys_mean"] <- mean(multi_sim_results$results[[i]]$esd$te)
    all_results[i, "esd_qalys_sd"] <- sd(multi_sim_results$results[[i]]$esd$te)
    all_results[i, "surgery_qalys_mean"] <- mean(multi_sim_results$results[[i]]$surgery$te)
    all_results[i, "surgery_qalys_sd"] <- sd(multi_sim_results$results[[i]]$surgery$te)
    all_results[i, "esd_costs_mean"] <- mean(multi_sim_results$results[[i]]$esd$tc)
    all_results[i, "esd_costs_sd"] <- sd(multi_sim_results$results[[i]]$esd$tc)
    all_results[i, "surgery_costs_mean"] <- mean(multi_sim_results$results[[i]]$surgery$tc)
    all_results[i, "surgery_costs_sd"] <- sd(multi_sim_results$results[[i]]$surgery$tc)
    all_results[i, "delta_qalys"] <- multi_sim_results$results[[i]]$delta_e
    all_results[i, "delta_costs"] <- multi_sim_results$results[[i]]$delta_c
    all_results[i, "icer"] <- multi_sim_results$results[[i]]$icer
  }
}

# Save the detailed results
write.csv(all_results, file = file.path(results_dir, "detailed_simulation_results.csv"), row.names = FALSE)

# Optional: Create a histogram of ICERs
if (require(ggplot2)) {
  icers <- sapply(multi_sim_results$results[1:multi_sim_results$completed_sims], 
                 function(x) if(x$completed) x$icer else NA)
  icers <- icers[!is.na(icers)]
  
  # Create histogram
  hist_data <- data.frame(icer = icers)
  p <- ggplot(hist_data, aes(x = icer)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = "Distribution of ICERs",
         x = "ICER (Cost per QALY)",
         y = "Frequency") +
    theme_minimal() +
    geom_vline(xintercept = mean(icers), color = "red", linetype = "dashed") +
    geom_vline(xintercept = median(icers), color = "blue", linetype = "dashed")
  
  # Save plot to file
  ggsave(file.path(results_dir, "icer_histogram.png"), p, width = 8, height = 6)
  
  # Create a scatterplot of costs vs QALYs
  if (nrow(all_results) > 0) {
    scatter_data <- data.frame(
      delta_qalys = all_results$delta_qalys,
      delta_costs = all_results$delta_costs
    )
    
    p2 <- ggplot(scatter_data, aes(x = delta_qalys, y = delta_costs)) +
      geom_point(color = "blue", alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      labs(title = "Cost-Effectiveness Plane",
           x = "Incremental QALYs",
           y = "Incremental Costs") +
      theme_minimal()
    
    # Add lines for different WTP thresholds
    wtp_thresholds <- c(20000, 50000, 100000)
    for (wtp in wtp_thresholds) {
      p2 <- p2 + geom_abline(slope = wtp, intercept = 0, 
                           linetype = "dotted", color = "red", alpha = 0.5)
    }
    
    # Save the cost-effectiveness plane
    ggsave(file.path(results_dir, "cost_effectiveness_plane.png"), p2, width = 8, height = 6)
    
    # Create boxplots for QALYs by strategy
    qalys_long <- data.frame(
      Strategy = rep(c("ESD", "Surgery"), each = nrow(all_results)),
      QALYs = c(all_results$esd_qalys_mean, all_results$surgery_qalys_mean)
    )
    
    p3 <- ggplot(qalys_long, aes(x = Strategy, y = QALYs, fill = Strategy)) +
      geom_boxplot() +
      labs(title = "QALYs by Strategy", 
           y = "QALYs",
           x = "") +
      theme_minimal()
    
    ggsave(file.path(results_dir, "qalys_boxplot.png"), p3, width = 8, height = 6)
  }
}

# Optional: Create a real-time monitoring function to track progress during long simulations
monitor_simulation_progress <- function(results_file = "./results/simulation_results_live.csv") {
  if (file.exists(results_file)) {
    # Read the current results
    results <- read.csv(results_file)
    
    # Print summary statistics
    cat("\nCurrent simulation progress:\n")
    cat("Completed simulations:", nrow(results), "\n")
    
    if (nrow(results) > 0) {
      cat("Mean ESD QALYs:", round(mean(results$esd_qalys), 3), "(SD:", round(sd(results$esd_qalys), 3), ")\n")
      cat("Mean Surgery QALYs:", round(mean(results$surgery_qalys), 3), "(SD:", round(sd(results$surgery_qalys), 3), ")\n")
      cat("Mean ESD Cost:", round(mean(results$esd_cost), 0), "(SD:", round(sd(results$esd_cost), 0), ")\n")
      cat("Mean Surgery Cost:", round(mean(results$surgery_cost), 0), "(SD:", round(sd(results$surgery_cost), 0), ")\n")
      cat("Mean ICER:", round(mean(results$icer), 0), "(SD:", round(sd(results$icer), 0), ")\n")
    }
  } else {
    cat("No simulation results file found at", results_file, "\n")
  }
}

# Make monitoring function available in global environment
assign("monitor_simulation_progress", monitor_simulation_progress, envir = .GlobalEnv)