# Fix source statement to use quotes
source("c:/Users/rma86/Desktop/gastro/BetaParmsFromQuantiles.R")
source("c:/Users/rma86/Desktop/gastro/GammaParmsFromQuantiles.R")

##################################### Functions ###########################################

# The MicroSim function for the simple microsimulation of the 'Sick-Sicker' model keeps track of what happens to each individual during each cycle. 

MicroSim <- function(v.M, n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, seed = 1, verbose = TRUE) {
  # Arguments:  
  # v.M:     vector of initial states for individuals 
  # n.i:     number of simulated individuals
  # n.t:     total number of cycles to run the model
  # v.n:     vector of health state names
  # d.c:     discount rate for costs
  # d.e:     discount rate for health outcome (QALYs)
  # TR.out:  should the output include a microsimulation trace? (default is TRUE)
  # TS.out:  should the output include a matrix of transitions between states? (default is TRUE)
  # seed:    starting seed number for random number generator (default is 1)
  # verbose: print detailed simulation progress information (default is TRUE)
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  
  v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e
  
  # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  # Initialize matrices to track time spent in S1 and S2 states
  m.S1_time <- m.S2_time <- matrix(0, nrow = n.i, ncol = n.t + 1,
                                   dimnames = list(paste("ind", 1:n.i, sep = " "),
                                                  paste("cycle", 0:n.t, sep = " ")))
  
  # Create transition frequency matrix to track how patients move between states
  m.trans_count <- matrix(0, nrow = n.s, ncol = n.s, 
                           dimnames = list(v.n, v.n))
  
  m.M[, 1] <- v.M                     # indicate the initial health state   
  
  # Print initial state distribution if verbose
  if (verbose) {
    cat("\nStarting simulation with", n.i, "individuals over", n.t, "cycles\n")
    cat("Initial state distribution:\n")
    print(table(factor(m.M[, 1], levels = v.n)) / n.i)
    cat("\n")
  }
  
  # Create matrices to track S1 and S2 substate distributions
  if (verbose) {
    s1_substate_tracker <- matrix(0, nrow = n.t+1, ncol = 9)
    s2_substate_tracker <- matrix(0, nrow = n.t+1, ncol = 9)
  }
    for (i in 1:n.i) {
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m.C[i, 1] <- Costs(m.M[i, 1], 0, 0)  # estimate costs per individual for the initial health state 
    m.E[i, 1] <- Effs(m.M[i, 1], 0, 0)  # estimate QALYs per individual for the initial health state
    
    for (t in 1:n.t) {
      # Print detailed cycle information for the first few individuals if verbose
      if (verbose && i <= 5 && t <= 5) {
        cat("Ind", i, "Cycle", t, ": Current state =", m.M[i, t], 
            "S1 time =", m.S1_time[i, t], "S2 time =", m.S2_time[i, t], "\n")
      }
      
      v.p <- Probs(m.M[i, t], m.S1_time[i, t], m.S2_time[i, t])  # calculate the transition probabilities at cycle t 
      
      # Print transition probabilities for the first few individuals if verbose
      if (verbose && i <= 3 && t <= 3) {
        cat("  Transition probabilities: ")
        for (s in 1:n.s) {
          if (v.p[s] > 0) {
            cat(v.n[s], "=", round(v.p[s], 4), " ")
          }
        }
        cat("\n")
      }
      
      m.M[i, t + 1] <- sample(v.n, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      
      # Track transitions between states
      from_state <- match(m.M[i, t], v.n)
      to_state <- match(m.M[i, t + 1], v.n)
      m.trans_count[from_state, to_state] <- m.trans_count[from_state, to_state] + 1
        m.C[i, t + 1] <- Costs(m.M[i, t + 1], m.S1_time[i, t + 1], m.S2_time[i, t + 1])   # estimate costs per individual during cycle t + 1
      m.E[i, t + 1] <- Effs(m.M[i, t + 1], m.S1_time[i, t + 1], m.S2_time[i, t + 1])   # estimate QALYs per individual during cycle t + 1
      
      # Update time counters for surveillance states
      if (m.M[i, t + 1] == "S1") {
        if (m.M[i, t] == "S1") {
          # If staying in S1, increment counter
          m.S1_time[i, t + 1] <- min(m.S1_time[i, t] + 1, 8)  # Cap at 8
        } else {
          # If newly entering S1, initialize counter to 1
          m.S1_time[i, t + 1] <- 1
        }
      }
      
      if (m.M[i, t + 1] == "S2") {
        if (m.M[i, t] == "S2") {
          # If staying in S2, increment counter
          m.S2_time[i, t + 1] <- min(m.S2_time[i, t] + 1, 8)  # Cap at 8
        } else {
          # If newly entering S2, initialize counter to 1
          m.S2_time[i, t + 1] <- 1
        }
      }
      
      # Print state transition for the first few individuals if verbose
      if (verbose && i <= 5 && t <= 5) {
        cat("  Transition:", m.M[i, t], "->", m.M[i, t + 1], "\n")
      }
      
    } # close the loop for the time points 
    
    # Progress indicator
    if(i/100 == round(i/100,0)) {          
      cat('\r', paste(i/n.i * 100, "% done", sep = " "))
    }
    
    # Print detailed information after every 1000 individuals if verbose
    if (verbose && i %% 1000 == 0) {
      cat("\nAfter simulating", i, "individuals:\n")
      
      # Track S1 and S2 substates after each individual
      s1_counts <- numeric(9)
      s2_counts <- numeric(9)
      
      for (s1_idx in 1:9) {
        s1_counts[s1_idx] <- sum(m.M[1:i, n.t + 1] == "S1" & m.S1_time[1:i, n.t + 1] == s1_idx)
      }
      
      for (s2_idx in 1:9) {
        s2_counts[s2_idx] <- sum(m.M[1:i, n.t + 1] == "S2" & m.S2_time[1:i, n.t + 1] == s2_idx)
      }
      
      # Print distribution of final health states
      cat("Current health state distribution at final cycle:\n")
      print(table(factor(m.M[1:i, n.t + 1], levels = v.n)) / i)
      cat("\nS1 substate distribution at final cycle:", s1_counts, "\n")
      cat("S2 substate distribution at final cycle:", s2_counts, "\n")
      
      # Calculate transition probabilities from counts
      cat("\nEmpirical transition matrix (proportion of transitions):\n")
      m.trans_prob <- m.trans_count / rowSums(m.trans_count)
      print(round(m.trans_prob, 4))
      cat("\n")
    }
    
  } # close the loop for the individuals 
  
  tc <- m.C %*% v.dwc       # total (discounted) cost per individual
  te <- m.E %*% v.dwe       # total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)        # average (discounted) cost 
  te_hat <- mean(te)        # average (discounted) QALYs
  
  if (TS.out == TRUE) {  # create a  matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns 
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) { # create a trace from the individual trajectories
    TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
    TR <- TR / n.i                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
    colnames(TR) <- v.n                                  # name the columns 
  } else {
    TR <- NULL
  }
  
  # Calculate S1 and S2 distributions by substate for final output
  if (verbose) {
    cat("\n\nFINAL SIMULATION RESULTS\n")
    cat("=======================\n\n")
    
    # Calculate distribution of S1 and S2 substates at the final cycle
    s1_substate_dist <- numeric(9)
    s2_substate_dist <- numeric(9)
    
    for (s1_idx in 1:9) {
      s1_substate_dist[s1_idx] <- sum(m.M[, n.t + 1] == "S1" & m.S1_time[, n.t + 1] == s1_idx)
    }
    
    for (s2_idx in 1:9) {
      s2_substate_dist[s2_idx] <- sum(m.M[, n.t + 1] == "S2" & m.S2_time[, n.t + 1] == s2_idx)
    }
    
    # Print final trace
    cat("State trace over time (proportion in each state by cycle):\n")
    print(round(TR, 4))
    
    # Print final state distribution
    cat("\nFinal health state distribution:\n")
    print(table(factor(m.M[, n.t + 1], levels = v.n)) / n.i)
    
    # Print final S1 and S2 substate distributions
    cat("\nFinal S1 substate distribution:", s1_substate_dist, "\n")
    cat("Final S2 substate distribution:", s2_substate_dist, "\n")
    
    # Print transition frequencies
    cat("\nTotal transitions between states:\n")
    print(m.trans_count)
      # Print empirical transition probabilities
    cat("\nEmpirical transition matrix (proportion of transitions):\n")
    # Handle division by zero (states with no transitions)
    row_sums <- rowSums(m.trans_count)
    m.trans_prob <- matrix(0, nrow = n.s, ncol = n.s, dimnames = list(v.n, v.n))
    for (i in 1:n.s) {
      if (row_sums[i] > 0) {
        m.trans_prob[i, ] <- m.trans_count[i, ] / row_sums[i]
      } else {
        # For rows with no transitions, indicate this with zeros instead of NaN
        m.trans_prob[i, ] <- 0
      }
    }
    print(round(m.trans_prob, 4))
  }
  
  # Store additional diagnostic information
  diag <- list(
    trans_count = m.trans_count,
    s1_time = m.S1_time,
    s2_time = m.S2_time
  )
  
  results <- list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, 
                  TS = TS, TR = TR, diag = diag) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it, S1_time_it = 0, S2_time_it = 0) { 
  # M_it:       health state occupied by individual i at cycle t (character variable)
  # S1_time_it: time spent in S1 state by individual i at cycle t (defaults to 0)
  # S2_time_it: time spent in S2 state by individual i at cycle t (defaults to 0)
  
  v.p.it <- rep(0, n.s)     # create vector of state transition probabilities with zeros
  names(v.p.it) <- v.n       # name the vector
  
  # Get index for time-dependent probabilities (capped at 8)
  S1_idx <- min(max(S1_time_it, 1), 8)  # Ensure index is between 1 and 8
  S2_idx <- min(max(S2_time_it, 1), 8)  # Ensure index is between 1 and 8
  
  # update v.p.it with the appropriate probabilities   
  if(M_it == "H1") {
    v.p.it["H1"] <- 0
    v.p.it["H2"] <- p.H1H2
    v.p.it["S1"] <- p.H1S1
    v.p.it["S2"] <- 0
    v.p.it["P"] <- 0
    v.p.it["D"] <- p.H1D
  } else if(M_it == "H2") {
    v.p.it["H1"] <- 0
    v.p.it["H2"] <- 0
    v.p.it["S1"] <- 0
    v.p.it["S2"] <- p.H2S2
    v.p.it["P"] <- p.H2P
    v.p.it["D"] <- p.H2D
  } else if(M_it == "S1") {
    v.p.it["H1"] <- 0
    v.p.it["H2"] <- p.S1H2[S1_idx]
    v.p.it["S1"] <- p.S1S1[S1_idx]
    v.p.it["S2"] <- 0
    v.p.it["P"] <- p.S1P[S1_idx]
    v.p.it["D"] <- p.S1D[S1_idx]
  } else if(M_it == "S2") {
    v.p.it["H1"] <- 0
    v.p.it["H2"] <- 0
    v.p.it["S1"] <- 0
    v.p.it["S2"] <- p.S2S2[S2_idx]
    v.p.it["P"] <- p.S2P[S2_idx]
    v.p.it["D"] <- p.S2D[S2_idx]
  } else if(M_it == "P") {
    v.p.it["H1"] <- 0
    v.p.it["H2"] <- 0
    v.p.it["S1"] <- 0
    v.p.it["S2"] <- 0
    v.p.it["P"] <- p.PP
    v.p.it["D"] <- p.PD
  } else if(M_it == "D") {
    v.p.it["H1"] <- 0
    v.p.it["H2"] <- 0
    v.p.it["S1"] <- 0
    v.p.it["S2"] <- 0
    v.p.it["P"] <- 0
    v.p.it["D"] <- 1
  }
  
  # Check if probabilities sum to 1 (allowing for small numerical errors)
  sum_p <- sum(v.p.it)
  if (abs(sum_p - 1) > 1e-6) {
    print(paste("Warning: Probabilities do not sum to 1 for state", M_it, 
                "S1_time =", S1_time_it, "S2_time =", S2_time_it, 
                "Sum =", sum_p))
    # Normalize probabilities to sum to 1
    v.p.it <- v.p.it / sum_p
  }
  
  return(v.p.it)
}


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function(M_it, S1_time_it = 0, S2_time_it = 0) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # S1_time_it: time spent in S1 state by individual i at cycle t (defaults to 0)
  # S2_time_it: time spent in S2 state by individual i at cycle t (defaults to 0)
  
  c.it <- 0                                  # by default the cost for everyone is zero 
  c.it[M_it == "H1"]  <- c.H1                # update the cost if H1
  c.it[M_it == "H2"]  <- c.H2                # update the cost if H2
  
  # For surveillance states, use the appropriate vector index based on time spent
  if (any(M_it == "S1")) {
    # Get index for time-dependent costs (capped at length of vector)
    S1_idx <- pmin(pmax(S1_time_it[M_it == "S1"], 1), length(c.S1))
    c.it[M_it == "S1"] <- c.S1[S1_idx]
  }
  
  if (any(M_it == "S2")) {
    # Get index for time-dependent costs (capped at length of vector)
    S2_idx <- pmin(pmax(S2_time_it[M_it == "S2"], 1), length(c.S2))
    c.it[M_it == "S2"] <- c.S2[S2_idx]
  }
  
  c.it[M_it == "P"]   <- c.P                 # update the cost if P
  
  return(c.it)                              # return the costs
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function(M_it, S1_time_it = 0, S2_time_it = 0, cl = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # S1_time_it: time spent in S1 state by individual i at cycle t (defaults to 0)
  # S2_time_it: time spent in S2 state by individual i at cycle t (defaults to 0)
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[M_it == "H1"]  <- u.H1    # update the utility if H1
  u.it[M_it == "H2"]  <- u.H2    # update the utility if H2
  
  # For surveillance states, use the appropriate vector index based on time spent
  if (any(M_it == "S1")) {
    # Get index for time-dependent utilities (capped at length of vector)
    S1_idx <- pmin(pmax(S1_time_it[M_it == "S1"], 1), length(u.S1))
    u.it[M_it == "S1"] <- u.S1[S1_idx]
  }
  
  if (any(M_it == "S2")) {
    # Get index for time-dependent utilities (capped at length of vector)
    S2_idx <- pmin(pmax(S2_time_it[M_it == "S2"], 1), length(u.S2))
    u.it[M_it == "S2"] <- u.S2[S2_idx]
  }
  
  u.it[M_it == "P"]  <- u.P      # update the utility if P
  u.it[M_it == "D"]  <- u.D      # update the utility if D
  
  QALYs <- u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}

### Helper functions 
# Collect all variables whose names start with "p."
get_all_p_vars <- function(env = parent.frame()) {
  p_names <- ls(env, pattern = "^p\\.")
  p_list <- setNames(lapply(p_names, function(nm) get(nm, envir = env)), p_names)
  return(p_list)
}

# Collect all variables whose names start with "c."
get_all_c_vars <- function(env = parent.frame()) {
  c_names <- ls(env, pattern = "^c\\.")
  c_list <- setNames(lapply(c_names, function(nm) get(nm, envir = env)), c_names)
  return(c_list)
}

# Collect all variables whose names start with "u."
get_all_u_vars <- function(env = parent.frame()) {
  u_names <- ls(env, pattern = "^u\\.")
  u_list <- setNames(lapply(u_names, function(nm) get(nm, envir = env)), u_names)
  return(u_list)
}

generate_beta_params <- function(p) {
  # Ensure p is within (0,1) and adjust quantiles to stay in (0,1)
  q1 <- max(min(p * 0.8, 1 - 1e-8), 1e-8)
  q2 <- max(min(p * 1.2, 1 - 1e-8), 1e-8)
  if (q1 > q2) { tmp <- q1; q1 <- q2; q2 <- tmp }
  beta.parms.from.quantiles(c(q1, q2))
}

# Function to generate a named list of random beta draws for all p.* variables
generate_random_p_draws <- function() {
  p_vars <- get_all_p_vars()
  random_draws <- list()
  
  for (name in names(p_vars)) {
    p_var <- p_vars[[name]]
    
    # Check if it's a vector probability
    if (is.vector(p_var) && length(p_var) > 1) {
      # Process vector probabilities element by element
      result_vector <- numeric(length(p_var))
      for (i in 1:length(p_var)) {
        p <- p_var[i]
        if (!is.numeric(p) || p <= 0 || p >= 1) {
          # Skip invalid probabilities
          result_vector[i] <- p
        } else {
          # Generate beta parameters and random draw
          tryCatch({
            params <- generate_beta_params(p)
            a <- if (!is.null(params$a)) params$a else params[1]
            b <- if (!is.null(params$b)) params$b else params[2]
            result_vector[i] <- rbeta(1, a, b)
          }, error = function(e) {
            # If error, keep original value
            result_vector[i] <- p
          })
        }
      }
      random_draws[[name]] <- result_vector
    } else {
      # Handle single value
      if (!is.numeric(p_var) || p_var <= 0 || p_var >= 1) {
        # Keep invalid values unchanged
        random_draws[[name]] <- p_var
      } else {
        # Generate beta parameters and random draw
        tryCatch({
          params <- generate_beta_params(p_var)
          a <- if (!is.null(params$a)) params$a else params[1]
          b <- if (!is.null(params$b)) params$b else params[2]
          random_draws[[name]] <- rbeta(1, a, b)
        }, error = function(e) {
          # If error, keep original value
          random_draws[[name]] <- p_var
        })
      }
    }
  }
  
  return(random_draws)
}

# Function to generate new random probability values and update global environment
generate_new_random_probs <- function(restore_originals = FALSE) {
  # Print debug information about global environment
  cat("Debugging generate_new_random_probs():\n")
  
  # Check for p.* variables in global environment
  p_vars <- get_all_p_vars(.GlobalEnv)
  cat("Found", length(p_vars), "p.* variables in global environment\n")
  
  if (length(p_vars) == 0) {
    cat("ERROR: No p.* variables found in global environment.\n")
    cat("Available variables in .GlobalEnv:\n")
    print(head(ls(.GlobalEnv), 20))
    return(list())
  }
  
  # Show first few p variables
  cat("First few p.* variables found:\n")
  print(head(names(p_vars), 5))
  
  # Generate random draws
  cat("Generating random draws...\n")
  random_draws <- generate_random_p_draws()
  
  cat("Generated", length(random_draws), "random draws\n")
  
  if (length(random_draws) > 0) {
    # Save original values
    original_values <- p_vars
    
    # Update global probability variables with random draws
    for (name in names(random_draws)) {
      assign(name, random_draws[[name]], envir = .GlobalEnv)
      cat("Assigned", name, "=", random_draws[[name]], "\n")
    }
    
    # Only restore original values if explicitly requested
    if (restore_originals) {
      for (name in names(original_values)) {
        assign(name, original_values[[name]], envir = .GlobalEnv)
      }
      cat("Restored original values\n")
    }
  } else {
    cat("ERROR: No random values generated\n")
  }
  
  # Return the random values that were generated
  return(random_draws)
}

generate_gamma_params <- function(val) {
  if (exists("gamma.parms.from.quantiles")) {
    q1 <- max(min(val * 0.8, 1 - 1e-8), 1e-8)  # Lower bound at 80% of value
    q2 <- max(min(val * 1.2, 1 - 1e-8), 1e-8)  # Upper bound at 120% of value
    if (q1 > q2) { tmp <- q1; q1 <- q2; q2 <- tmp }
    return(gamma.parms.from.quantiles(c(q1, q2)))
  } else {
    warning("gamma.parms.from.quantiles function not found")
    return(NA)
  }
}

# Function to generate random draws for cost variables
generate_random_c_draws <- function() {
  c_vars <- get_all_c_vars()
  
  # Apply generate_gamma_params to each value
  gamma_params_list <- lapply(c_vars, generate_gamma_params)
  
  # Generate random gamma values using these parameters
  random_draws <- mapply(function(parms, name) {
    if(is.na(parms)[1]) {
      return(c_vars[[name]]) # Return original if NA
    }
    a <- if(!is.null(parms$shape)) parms$shape else parms[1]
    scale <- if(!is.null(parms$scale)) parms$scale else parms[2]
    rgamma(1, shape = a, scale = scale)
  }, gamma_params_list, names(gamma_params_list), SIMPLIFY = FALSE)
  
  return(random_draws)
}

# Function to generate random draws for utility variables
generate_random_u_draws <- function() {
  u_vars <- get_all_u_vars()
  
  # For utilities, we need to ensure values stay between 0 and 1
  # We'll use a beta distribution for utilities
  random_draws <- lapply(u_vars, function(u_val) {
    # Skip vectors, only process single values
    if(length(u_val) > 1) return(u_val)
    if(u_val <= 0 || u_val >= 1) return(u_val) # Return as is if out of (0,1) range
    
    # Use beta distribution centered at u_val with small variance
    alpha <- u_val * 10
    beta <- (1-u_val) * 10
    rbeta(1, alpha, beta)
  })
  
  return(random_draws)
}

# Function to update cost variables with random values
generate_new_random_c_values <- function(return_originals = FALSE) {
  # Generate random draws
  random_c_draws <- generate_random_c_draws()
  
  # Save original values
  original_values <- get_all_c_vars()
  
  # Update global variables
  for(name in names(random_c_draws)) {
    # Don't update vector variables
    if(length(original_values[[name]]) == 1) {
      assign(name, random_c_draws[[name]], envir = .GlobalEnv)
    }
  }
  
  # Return original values without printing to console if needed
  if(return_originals) {
    return(original_values)
  } else {
    # Return TRUE to indicate function completed successfully
    return(TRUE)
  }
}

# Function to update utility variables with random values
generate_new_random_u_values <- function(return_originals = FALSE) {
  # Generate random draws
  random_u_draws <- generate_random_u_draws()
  
  # Save original values
  original_values <- get_all_u_vars()
  
  # Update global variables
  for(name in names(random_u_draws)) {
    # Don't update vector variables
    if(length(original_values[[name]]) == 1) {
      assign(name, random_u_draws[[name]], envir = .GlobalEnv)
    }
  }
  
  # Return original values without printing to console if needed
  if(return_originals) {
    return(original_values)
  } else {
    # Return TRUE to indicate function completed successfully
    return(TRUE)
  }
}

# Function to initialize model parameters in the global environment
initialize_model_parameters <- function() {
  # Transition probabilities (per cycle)
  assign("p.H1D", 0.0050, envir = .GlobalEnv)              # 3 month mortality rate after ESD)
  assign("p.H1H2", 0.0886, envir = .GlobalEnv)              # probability of switching to surgery after ESD
  
  # Calculate and assign dependent probabilities
  assign("p.H1S1", 1 - 0.0050 - 0.0886, envir = .GlobalEnv)  # probability of staying in surveillance after ESD
  assign("p.H2D", 0.0463, envir = .GlobalEnv)             # short-term mortality rate after surgery
  assign("p.H2P", 0.0082, envir = .GlobalEnv)             # probability of clinical failure after ESD
  assign("p.H2S2", 1 - 0.0463 - 0.0082, envir = .GlobalEnv)  # probability of staying in surveillance after surgery
 
  # S1 state transition probabilities
  assign("p.S1H2", c(0, 0.03, 0.03, 0.0249, 0.0249, 0.0249, 0.0249, 0.02), envir = .GlobalEnv)  
  assign("p.S1P", rep(0, 8), envir = .GlobalEnv)     
  assign("p.S1D", c(0.002, 0.0034, 0.0034, 0.0028, 0.0028, 0.0028, 0.0028, 0.0023), envir = .GlobalEnv)
  
  # Calculate p.S1S1 for each time point
  p.S1S1_values <- 1 - get("p.S1H2", envir = .GlobalEnv) - get("p.S1P", envir = .GlobalEnv) - get("p.S1D", envir = .GlobalEnv)
  assign("p.S1S1", p.S1S1_values, envir = .GlobalEnv)
  
  # S2 state transition probabilities
  assign("p.S2P", c(0.0000, 0.0064, 0.0064, 0.0064, 0.0038, 0.0038, 0.0038, 0.0025), envir = .GlobalEnv)
  assign("p.S2D", c(0.0058, 0.0058, 0.0044, 0.0044, 0.0026, 0.0026, 0.0026, 0.0018), envir = .GlobalEnv)
  
  # Calculate p.S2S2 for each time point
  p.S2S2_values <- 1 - get("p.S2P", envir = .GlobalEnv) - get("p.S2D", envir = .GlobalEnv)
  assign("p.S2S2", p.S2S2_values, envir = .GlobalEnv)
  
  # Palliative care probabilities
  assign("p.PD", 0.5, envir = .GlobalEnv)  # assume patient dies with 50% chance when in palliative care
  assign("p.PP", 1 - 0.5, envir = .GlobalEnv)  # probability of staying in palliative care
  
  # Cost and utility inputs (we should initialize these too)
  assign("c.H1", 18347.00, envir = .GlobalEnv)                 # cost of ESD
  assign("c.H2", 80000.00, envir = .GlobalEnv)                # cost of surgery
  assign("c.S1", c(4137.50, 2068.75, 2068.75, 1443.75, 1443.75, 1443.75, 1443.75, 1443.75), envir = .GlobalEnv) 
  assign("c.S2", c(3275.00, 1637.50, 1637.50, 818.75, 818.75, 818.75, 818.75, 818.75), envir = .GlobalEnv)
  assign("c.P", 9600, envir = .GlobalEnv)                  # cost of palliative

  assign("u.H1", 0.7, envir = .GlobalEnv)                 # utility when PreESD 
  assign("u.H2", 0.7, envir = .GlobalEnv)                 # utility when PreSurgery 
  assign("u.S1", c(0.73, 0.75, 0.75, 0.76, 0.76, 0.76, 0.76, 0.76), envir = .GlobalEnv)
  assign("u.S2", c(0.73, 0.75, 0.75, 0.76, 0.76, 0.76, 0.76, 0.76), envir = .GlobalEnv)
  assign("u.P", 0.5, envir = .GlobalEnv)                 # utility when palliative
  assign("u.D", 0, envir = .GlobalEnv)                   # utility when dead
  
  # Model structure parameters  
  assign("n.s", 6, envir = .GlobalEnv)  # number of health states
  assign("v.n", c("H1","H2","S1","S2","P","D"), envir = .GlobalEnv)  # the model states
  
  cat("Model parameters initialized successfully.\n")
  cat("Probability (p.*) variables:", length(ls(pattern = "^p\\.", envir = .GlobalEnv)), "\n")
  cat("Cost (c.*) variables:", length(ls(pattern = "^c\\.", envir = .GlobalEnv)), "\n")
  cat("Utility (u.*) variables:", length(ls(pattern = "^u\\.", envir = .GlobalEnv)), "\n")
}