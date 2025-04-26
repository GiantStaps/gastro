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

##################################### Unified simulation function ##########################
run_simulation <- function(n.i = 1000, n.t = 20, sensitivity = FALSE, time_limit_seconds = NULL, 
                           verbose = FALSE, sim_index = 1) {
  # Start timing if time limit is specified
  start_time <- Sys.time()
  
  # Initialize model parameters
  initialize_model_parameters(force = TRUE)
  
  # Set seed for reproducibility
  set.seed(12345 + sim_index)
  
  # If running sensitivity analysis, generate random values for probabilities, costs, and utilities
  if (sensitivity) {
    generate_new_random_probs(print_output = FALSE)
    generate_new_random_c_values(print_output = FALSE)
    generate_new_random_u_values(print_output = FALSE)
  }

  # Define the initial state vectors
  v.M_1 <- rep("H1", n.i)  # all start in the H1 (ESD) state
  v.M_2 <- rep("H2", n.i)  # all start in the H2 (Surgery) state
  
  # Treatment names
  v.Trt <- c("ESD", "Surgery")
  
  # Discount rates (annualized)
  d.c <- d.e <- 0.03 / 4  # discount rate for costs and QALYs per quarter
  
  # Run ESD simulation
  if (verbose) cat("Running ESD simulation:\n")
  sim_esd <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = (12345 + sim_index*2))
  
  # Check if we've exceeded the time limit (if specified)
  if (!is.null(time_limit_seconds) && 
      difftime(Sys.time(), start_time, units = "secs") > time_limit_seconds) {
    if (verbose) cat("Time limit reached, returning partial results\n")
    # Divide QALYs by 4 for annualization
    sim_esd$te <- sim_esd$te / 4
    return(list(
      esd = sim_esd,
      surgery = NULL,
      completed = FALSE,
      runtime = difftime(Sys.time(), start_time, units = "secs")
    ))
  }
  
  # Run Surgery simulation
  if (verbose) cat("Running Surgery simulation:\n")
  sim_surgery <- MicroSim(v.M_2, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = (12345 + sim_index*2 + 1))
  
  # Divide QALYs by 4 for annualization
  sim_esd$te <- sim_esd$te / 4
  sim_surgery$te <- sim_surgery$te / 4
  
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
run_multiple_simulations <- function(n_reps = 100, sensitivity = FALSE, time_limit_seconds = NULL, 
                                   n.i = 1000, n.t = 20, verbose = FALSE, seed = 12345) {
  # Prepare storage for results
  results <- list()
  completed_sims <- 0
  total_start_time <- Sys.time()
  
  # Print simulation parameters
  simulation_type <- if(sensitivity) "sensitivity analysis" else "baseline"
  time_limit_msg <- if(!is.null(time_limit_seconds)) paste("with time limit of", time_limit_seconds, "seconds each") else "with no time limit"
  
  cat("Starting", n_reps, simulation_type, "simulations", time_limit_msg, "\n")
  cat("Total individuals per simulation:", n.i, "\n")
  cat("Total cycles per simulation:", n.t, "\n\n")
  
  # Create results directory if it doesn't exist
  results_dir <- "./results"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  
  # Run simulations
  for (i in 1:n_reps) {
    cat("Starting simulation", i, "of", n_reps, "\n")
    sim_start_time <- Sys.time()
    
    # Run a single simulation with appropriate parameters
    sim_result <- run_simulation(
      n.i = n.i, 
      n.t = n.t, 
      sensitivity = sensitivity,
      time_limit_seconds = time_limit_seconds, 
      verbose = verbose, 
      sim_index = seed + i - 1
    )
    
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
      # Calculate mean costs and QALYs for reporting (QALYs already annualized)
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
          
      # Add time taken to the result for UI display
      sim_result$time_taken <- as.numeric(sim_runtime)
      
      # Call the updateIncrementalResults function if it exists in parent frame
      if (exists("updateIncrementalResults", envir = parent.frame())) {
        updateIncrementalResults_fn <- get("updateIncrementalResults", envir = parent.frame())
        updateIncrementalResults_fn(i, sim_result)
      }
    } else {
      cat(" (incomplete - time limit reached)\n")
    }
  }
  
  total_runtime <- difftime(Sys.time(), total_start_time, units = "secs")
  cat("\nAll simulations completed in", round(total_runtime, 2), "seconds\n")
  cat("Completed simulations:", completed_sims, "out of", n_reps, "\n\n")
  
  # Compile summary statistics for completed simulations
  if (completed_sims > 0) {
    # Filter out incomplete simulations
    completed_results <- Filter(function(x) x$completed, results)
    
    # Extract metrics
    icers <- sapply(completed_results, function(x) x$icer)
    delta_costs <- sapply(completed_results, function(x) x$delta_c)
    delta_qalys <- sapply(completed_results, function(x) x$delta_e)
    esd_costs <- sapply(completed_results, function(x) mean(x$esd$tc))
    surgery_costs <- sapply(completed_results, function(x) mean(x$surgery$tc))
    esd_qalys <- sapply(completed_results, function(x) mean(x$esd$te))
    surgery_qalys <- sapply(completed_results, function(x) mean(x$surgery$te))
    
    # Print summary statistics
    cat("Summary of completed simulations:\n")
    cat("Mean ICER:", round(mean(icers), 2), "(SD:", round(sd(icers), 2), ")\n")
    cat("Mean Incremental Cost:", round(mean(delta_costs), 2), "(SD:", round(sd(delta_costs), 2), ")\n")
    cat("Mean Incremental QALYs:", round(mean(delta_qalys), 4), "(SD:", round(sd(delta_qalys), 4), ")\n")
    cat("Mean ESD Cost:", round(mean(esd_costs), 2), "(SD:", round(sd(esd_costs), 2), ")\n")
    cat("Mean Surgery Cost:", round(mean(surgery_costs), 2), "(SD:", round(sd(surgery_costs), 2), ")\n")
    cat("Mean ESD QALYs:", round(mean(esd_qalys), 4), "(SD:", round(sd(esd_qalys), 4), ")\n")
    cat("Mean Surgery QALYs:", round(mean(surgery_qalys), 4), "(SD:", round(sd(surgery_qalys), 4), ")\n")
    
    # Calculate percentiles for ICER if we have enough simulations
    if (length(icers) >= 4) {
      icer_percentiles <- quantile(icers, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      cat("\nICER Percentiles:\n")
      cat("2.5%:", round(icer_percentiles[1], 2), "\n")
      cat("25%:", round(icer_percentiles[2], 2), "\n")
      cat("50% (median):", round(icer_percentiles[3], 2), "\n")
      cat("75%:", round(icer_percentiles[4], 2), "\n")
      cat("97.5%:", round(icer_percentiles[5], 2), "\n")
    }
  }
  
  # Return a structured result for further processing
  return(list(
    results = results,
    completed_sims = completed_sims,
    total_runtime = total_runtime,
    simulation_type = if(sensitivity) "sensitivity" else "baseline",
    n_reps = n_reps,
    n_individuals = n.i,
    n_cycles = n.t,
    summary = if(completed_sims > 0) {
      list(
        mean_icer = mean(icers),
        sd_icer = sd(icers),
        mean_delta_cost = mean(delta_costs),
        sd_delta_cost = sd(delta_costs),
        mean_delta_qaly = mean(delta_qalys),
        sd_delta_qaly = sd(delta_qalys),
        mean_esd_cost = mean(esd_costs),
        sd_esd_cost = sd(esd_costs),
        mean_surgery_cost = mean(surgery_costs),
        sd_surgery_cost = sd(surgery_costs),
        mean_esd_qaly = mean(esd_qalys),
        sd_esd_qaly = sd(esd_qalys),
        mean_surgery_qaly = mean(surgery_qalys),
        sd_surgery_qaly = sd(surgery_qalys),
        icer_percentiles = if(length(icers) >= 4) quantile(icers, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) else NULL
      )
    } else NULL
  ))
}

# Convenience functions for backward compatibility
run_baseline_simulation <- function(n.i = 1000, n.t = 20, verbose = FALSE, n_reps = 1, seed = 12345) {
  # For single replication, return just that result
  if (n_reps == 1) {
    return(run_simulation(n.i, n.t, sensitivity = FALSE, verbose = verbose, sim_index = seed))
  } else {
    return(run_multiple_simulations(n_reps, sensitivity = FALSE, n.i = n.i, n.t = n.t, 
                                   verbose = verbose, seed = seed))
  }
}

sensitivity_analysis <- function(n_reps = 1000, time_limit_seconds = 60, n.i = 1000, n.t = 20, 
                                verbose = FALSE, seed = 12345) {
  return(run_multiple_simulations(n_reps, sensitivity = TRUE, time_limit_seconds = time_limit_seconds,
                                n.i = n.i, n.t = n.t, verbose = verbose, seed = seed))
}