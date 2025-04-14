# Reset environment and initialize
rm(list = ls())  # Clear environment first

# First source the functions.r file explicitly
source("c:/Users/rma86/Desktop/gastro/functions.r")  # Use lowercase .r to match actual file
source("c:/Users/rma86/Desktop/gastro/BetaParmsFromQuantiles.R")
source("c:/Users/rma86/Desktop/gastro/GammaParmsFromQuantiles.R")

# Load the function definitions from the microsimulation model
source("c:/Users/rma86/Desktop/gastro/microsimulation_model.R")

# Test just initialize_model_parameters
if (exists("initialize_model_parameters")) {
  cat("initialize_model_parameters function exists\n")
  initialize_model_parameters(force = TRUE)
} else {
  cat("ERROR: initialize_model_parameters function not found!\n")
}

# Test just one simulation - should run faster than the full set
if (exists("run_timed_simulation")) {
  cat("\nTesting a single simulation run...\n")
  single_sim <- run_timed_simulation(time_limit_seconds = 60, n.i = 100, n.t = 10, verbose = TRUE)
  if (!is.null(single_sim)) {
    cat("\nSingle simulation completed successfully!\n")
    if (single_sim$completed) {
      cat("ICER:", round(single_sim$icer, 2), "\n")
    }
  }
} else {
  cat("ERROR: run_timed_simulation function not found!\n")
}