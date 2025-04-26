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
run_timed_simulation <- function(time_limit_seconds = 60, n.i = 1000, n.t = 20, verbose = FALSE, sim_index = 1) {
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
  
  # Discount rates (annualized)
  d.c <- d.e <- 0.03 / 4  # discount rate for costs and QALYs per quarter
  
  # Run ESD and Surgery simulations
  if (verbose) cat("Running ESD simulation:\n")
  sim_esd <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = (12345 + sim_index*2))
  
  # Check if we've exceeded the time limit
  if (difftime(Sys.time(), start_time, units = "secs") > time_limit_seconds) {
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
sensitivity_analysis <- function(n_reps = 1000, time_limit_seconds = 60, n.i = 1000, n.t = 20, verbose = FALSE, seed = 12345) {
  # Prepare storage for results
  results <- list()
  completed_sims <- 0
  total_start_time <- Sys.time()
  
  cat("Starting", n_reps, "simulations with time limit of", time_limit_seconds, "seconds each\n")
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
  
  for (i in 1:n_reps) {
    cat("Starting simulation", i, "of", n_reps, "\n")
    sim_start_time <- Sys.time()
    
    # Run a single timed simulation with consistent seeding
    sim_result <- run_timed_simulation(time_limit_seconds, n.i, n.t, verbose, sim_index = seed + i - 1)
    
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
  cat("Completed simulations:", completed_sims, "out of", n_reps, "\n\n")
  
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

# New: Run baseline simulation with initialized model input (no randomization)
run_baseline_simulation <- function(n.i = 1000, n.t = 20, verbose = FALSE, n_reps = 1, seed = 12345) {
  results <- vector("list", n_reps)
  for (i in 1:n_reps) {
    # Set seed at the beginning of each replication
    set.seed(seed + i - 1)
    initialize_model_parameters(force = TRUE)
    v.M_1 <- rep("H1", n.i)
    v.M_2 <- rep("H2", n.i)
    d.c <- d.e <- 0.03 / 4
    # Pass the seed to MicroSim but don't set it there
    sim_esd <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = seed + i*2)
    sim_surgery <- MicroSim(v.M_2, n.i, n.t, v.n, d.c, d.e, verbose = verbose, seed = seed + i*2 + 1)
    sim_esd$te <- sim_esd$te / 4
    sim_surgery$te <- sim_surgery$te / 4
    delta_c <- mean(sim_surgery$tc) - mean(sim_esd$tc)
    delta_e <- mean(sim_surgery$te) - mean(sim_esd$te)
    icer <- delta_c / delta_e
    results[[i]] <- list(
      esd = sim_esd,
      surgery = sim_surgery,
      delta_c = delta_c,
      delta_e = delta_e,
      icer = icer
    )
  }
  if (n_reps == 1) {
    return(results[[1]])
  } else {
    esd_qalys <- sapply(results, function(res) mean(res$esd$te))
    surgery_qalys <- sapply(results, function(res) mean(res$surgery$te))
    esd_costs <- sapply(results, function(res) mean(res$esd$tc))
    surgery_costs <- sapply(results, function(res) mean(res$surgery$tc))
    delta_qalys <- sapply(results, function(res) res$delta_e)
    delta_costs <- sapply(results, function(res) res$delta_c)
    icers <- sapply(results, function(res) res$icer)
    return(list(
      n_reps = n_reps,
      esd_qalys_mean = mean(esd_qalys), esd_qalys_sd = sd(esd_qalys),
      surgery_qalys_mean = mean(surgery_qalys), surgery_qalys_sd = sd(surgery_qalys),
      esd_costs_mean = mean(esd_costs), esd_costs_sd = sd(esd_costs),
      surgery_costs_mean = mean(surgery_costs), surgery_costs_sd = sd(surgery_costs),
      delta_qalys_mean = mean(delta_qalys), delta_qalys_sd = sd(delta_qalys),
      delta_costs_mean = mean(delta_costs), delta_costs_sd = sd(delta_costs),
      icer_mean = mean(icers), icer_sd = sd(icers),
      all_results = results
    ))
  }
}