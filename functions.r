# Fix source statement to use quotes
source("BetaParmsFromQuantiles.R")
source("GammaParmsFromQuantiles.R")

##################################### Functions ###########################################

MicroSim <- function(v.M, n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, seed = 1, verbose = TRUE, output_json = NULL) {
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
    
    # Progress indicator - MODIFIED TO ONLY SHOW IN VERBOSE MODE
    if(verbose && i/100 == round(i/100,0)) {          
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
  
  # Per-cycle logging (after all individuals for each cycle)
  if (!is.null(output_json)) {
    for (cycle in 1:(n.t + 1)) {
      # State distribution
      state_counts <- as.list(table(factor(m.M[, cycle], levels = v.n)))
      names(state_counts) <- v.n
      # Costs and QALYs up to this cycle
      costs <- rowSums(m.C[, 1:cycle, drop = FALSE])
      qalys <- rowSums(m.E[, 1:cycle, drop = FALSE])
      # Statistics
      cost_avg <- mean(costs)
      cost_sd <- sd(costs)
      qaly_avg <- mean(qalys)
      qaly_sd <- sd(qalys)
      # ICER: not meaningful for a single arm, so set to NA here
      icer_avg <- NA
      icer_sd <- NA
      # Prepare JSON object
      cycle_result <- list(
        state_counts = state_counts,
        cost = list(avg = cost_avg, std = cost_sd),
        qaly = list(avg = qaly_avg, std = qaly_sd),
        icer = list(avg = icer_avg, std = icer_sd)
      )
      # Append to file as JSON line
      json_line <- jsonlite::toJSON(cycle_result, auto_unbox = TRUE)
      write(json_line, file = output_json, append = TRUE)
    }
  }
  
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
get_all_p_vars <- function(env = .GlobalEnv) {  # changed default env from parent.frame() to .GlobalEnv
  p_names <- ls(env, pattern = "^p\\.")
  p_list <- setNames(lapply(p_names, function(nm) get(nm, envir = env)), p_names)
  return(p_list)
}

# Fix get_all_c_vars to use .GlobalEnv by default like get_all_p_vars
get_all_c_vars <- function(env = .GlobalEnv) {
  c_names <- ls(env, pattern = "^c\\.")
  c_list <- setNames(lapply(c_names, function(nm) get(nm, envir = env)), c_names)
  return(c_list)
}

# Collect all variables whose names start with "u."
get_all_u_vars <- function(env = .GlobalEnv) {
  u_names <- ls(env, pattern = "^u\\.")
  u_list <- setNames(lapply(u_names, function(nm) get(nm, envir = env)), u_names)
  return(u_list)
}

# Revised generate_beta_params: cap q2 at 0.999 to avoid numerical issues
generate_beta_params <- function(p) {
  # Ensure p is strictly in (0,1)
  if(p <= 0) p <- 1e-8
  if(p >= 1) p <- 1 - 1e-8
  
  # Calculate quantiles, but cap q2 at 0.999
  q1 <- max(min(p * 0.8, 1 - 1e-8), 1e-8)
  q2 <- max(min(p * 1.2, 0.999), 1e-8)  # Cap at 0.999 instead of 1-1e-8
  
  if(q1 > q2) { tmp <- q1; q1 <- q2; q2 <- tmp }
  beta.parms.from.quantiles(c(q1, q2))
}

# Modified function to generate a named list of random beta draws for all p.* variables
# with constraints for dependent probabilities
generate_random_p_draws <- function(print_output = FALSE) {
  p_vars <- get_all_p_vars()
  if(print_output) cat("Processing", length(p_vars), "probability variables\n")
  random_draws <- list()
  
  # Define constrained variables that should not be randomly generated
  constrained_vars <- c("p.H1S1", "p.H2S2", "p.S1S1", "p.S2S2", "p.PP")
  
  for (name in names(p_vars)) {
    # Skip constrained variables - they will be calculated later
    if (name %in% constrained_vars) {
      if(print_output) cat("Skipping constrained variable:", name, "\n")
      next
    }
    
    p_var <- p_vars[[name]]
    
    # Print variable information only if requested
    if(print_output) {
      if (is.vector(p_var) && length(p_var) > 1) {
        cat("Variable:", name, "is a vector with", length(p_var), "elements:", head(p_var), "...\n")
      } else {
        cat("Variable:", name, "=", p_var, "\n")
      }
    }
    
    # Handle vector probabilities (such as p.S1H2, p.S1P, p.S1D)
    if (is.vector(p_var) && length(p_var) > 1) {
      # Process vector probabilities element by element
      result_vector <- numeric(length(p_var))
      for (i in 1:length(p_var)) {
        p <- p_var[i]
        # Explicitly check for p = 0
        if (!is.numeric(p) || p == 0 || p >= 1) {
          if(print_output) cat("  Element", i, "value", p, "is zero or invalid, keeping original\n")
          result_vector[i] <- p
        } else {
          # Generate beta parameters and random draw
          tryCatch({
            params <- generate_beta_params(p)
            a <- if (!is.null(params$a)) params$a else params[1]
            b <- if (!is.null(params$b)) params$b else params[2]
            result_vector[i] <- rbeta(1, a, b)
            if(print_output) cat("  Element", i, "value", p, "-> random draw:", result_vector[i], "\n")
          }, error = function(e) {
            # Always print errors
            cat("  Error for element", i, ":", e$message, "\n")
            result_vector[i] <- p
          })
        }
      }
      random_draws[[name]] <- result_vector
      if(print_output) cat("  Saved vector result for", name, "\n")
    } else {
      # Handle single value
      # Explicitly check for p = 0
      if (!is.numeric(p_var) || p_var == 0 || p_var >= 1) {
        if(print_output) cat("  Value", p_var, "is zero or invalid, keeping original\n")
        random_draws[[name]] <- p_var
      } else {
        # Generate beta parameters and random draw
        tryCatch({
          params <- generate_beta_params(p_var)
          a <- if (!is.null(params$a)) params$a else params[1]
          b <- if (!is.null(params$b)) params$b else params[2]
          draw <- rbeta(1, a, b)
          random_draws[[name]] <- draw
          if(print_output) {
            cat("  Generated parameters: a =", a, "b =", b, "\n")
            cat("  Random draw:", draw, "\n")
          }
        }, error = function(e) {
          # Always print errors
          cat("  Error:", e$message, "\n")
          random_draws[[name]] <- p_var
        })
      }
    }
  }
  
  # Now calculate constrained variables based on their dependencies
  
  # 1. Calculate p.H1S1 = 1 - p.H1D - p.H1H2
  if ("p.H1D" %in% names(random_draws) && "p.H1H2" %in% names(random_draws)) {
    p_h1d <- random_draws[["p.H1D"]]
    p_h1h2 <- random_draws[["p.H1H2"]]
    
    # Ensure the sum of probabilities doesn't exceed 1
    if (p_h1d + p_h1h2 > 0.999) {
      # Scale down both probabilities to make room for p.H1S1
      scale_factor <- 0.9 / (p_h1d + p_h1h2)
      p_h1d <- p_h1d * scale_factor
      p_h1h2 <- p_h1h2 * scale_factor
      random_draws[["p.H1D"]] <- p_h1d
      random_draws[["p.H1H2"]] <- p_h1h2
      if(print_output) cat("Rescaled p.H1D and p.H1H2 to ensure sum < 1\n")
    }
    
    random_draws[["p.H1S1"]] <- 1 - p_h1d - p_h1h2
    if(print_output) cat("Calculated constrained p.H1S1 =", random_draws[["p.H1S1"]], "\n")
  }
  
  # 2. Calculate p.H2S2 = 1 - p.H2P - p.H2D
  if ("p.H2P" %in% names(random_draws) && "p.H2D" %in% names(random_draws)) {
    p_h2p <- random_draws[["p.H2P"]]
    p_h2d <- random_draws[["p.H2D"]]
    
    # Ensure the sum of probabilities doesn't exceed 1
    if (p_h2p + p_h2d > 0.999) {
      # Scale down both probabilities to make room for p.H2S2
      scale_factor <- 0.9 / (p_h2p + p_h2d)
      p_h2p <- p_h2p * scale_factor
      p_h2d <- p_h2d * scale_factor
      random_draws[["p.H2P"]] <- p_h2p
      random_draws[["p.H2D"]] <- p_h2d
      if(print_output) cat("Rescaled p.H2P and p.H2D to ensure sum < 1\n")
    }
    
    random_draws[["p.H2S2"]] <- 1 - p_h2p - p_h2d
    if(print_output) cat("Calculated constrained p.H2S2 =", random_draws[["p.H2S2"]], "\n")
  }
  
  # 3. Calculate p.S1S1 = 1 - p.S1H2 - p.S1P - p.S1D (vector)
  if (all(c("p.S1H2", "p.S1P", "p.S1D") %in% names(random_draws))) {
    p_s1h2 <- random_draws[["p.S1H2"]]
    p_s1p <- random_draws[["p.S1P"]]
    p_s1d <- random_draws[["p.S1D"]]
    
    # For vector probabilities, calculate p.S1S1 element-wise
    p_s1s1 <- numeric(length(p_s1h2))
    
    for (i in 1:length(p_s1h2)) {
      # Ensure the sum of vector probabilities doesn't exceed 1 for any element
      sum_probs <- p_s1h2[i] + p_s1p[i] + p_s1d[i]
      if (sum_probs > 0.999) {
        # Scale down probabilities to make room for p.S1S1
        scale_factor <- 0.9 / sum_probs
        p_s1h2[i] <- p_s1h2[i] * scale_factor
        p_s1p[i] <- p_s1p[i] * scale_factor
        p_s1d[i] <- p_s1d[i] * scale_factor
        if(print_output) cat("Rescaled p.S1H2[", i, "], p.S1P[", i, "], p.S1D[", i, "] to ensure sum < 1\n", sep="")
      }
      p_s1s1[i] <- 1 - p_s1h2[i] - p_s1p[i] - p_s1d[i]
    }
    
    # Update the random draws with potentially rescaled values
    random_draws[["p.S1H2"]] <- p_s1h2
    random_draws[["p.S1P"]] <- p_s1p
    random_draws[["p.S1D"]] <- p_s1d
    random_draws[["p.S1S1"]] <- p_s1s1
    
    if(print_output) cat("Calculated constrained p.S1S1 vector\n")
  }
  
  # 4. Calculate p.S2S2 = 1 - p.S2P - p.S2D (vector)
  if (all(c("p.S2P", "p.S2D") %in% names(random_draws))) {
    p_s2p <- random_draws[["p.S2P"]]
    p_s2d <- random_draws[["p.S2D"]]
    
    # For vector probabilities, calculate p.S2S2 element-wise
    p_s2s2 <- numeric(length(p_s2p))
    
    for (i in 1:length(p_s2p)) {
      # Ensure the sum of vector probabilities doesn't exceed 1 for any element
      sum_probs <- p_s2p[i] + p_s2d[i]
      if (sum_probs > 0.999) {
        # Scale down probabilities to make room for p.S2S2
        scale_factor <- 0.9 / sum_probs
        p_s2p[i] <- p_s2p[i] * scale_factor
        p_s2d[i] <- p_s2d[i] * scale_factor
        if(print_output) cat("Rescaled p.S2P[", i, "], p.S2D[", i, "] to ensure sum < 1\n", sep="")
      }
      p_s2s2[i] <- 1 - p_s2p[i] - p_s2d[i]
    }
    
    # Update the random draws with potentially rescaled values
    random_draws[["p.S2P"]] <- p_s2p
    random_draws[["p.S2D"]] <- p_s2d
    random_draws[["p.S2S2"]] <- p_s2s2
    
    if(print_output) cat("Calculated constrained p.S2S2 vector\n")
  }
  
  # 5. Calculate p.PP = 1 - p.PD
  if ("p.PD" %in% names(random_draws)) {
    p_pd <- random_draws[["p.PD"]]
    
    # Ensure p.PD doesn't exceed 1
    if (p_pd > 0.999) {
      p_pd <- 0.999
      random_draws[["p.PD"]] <- p_pd
      if(print_output) cat("Capped p.PD at 0.999 to ensure p.PP ≥ 0.001\n")
    }
    
    random_draws[["p.PP"]] <- 1 - p_pd
    if(print_output) cat("Calculated constrained p.PP =", random_draws[["p.PP"]], "\n")
  }
  
  if(print_output) cat("Final random_draws list contains", length(random_draws), "items\n")
  return(random_draws)
}

# Function to generate new random probability values and update global environment
generate_new_random_probs <- function(restore_originals = FALSE, print_output = FALSE) {
  # Check for p.* variables in global environment
  p_vars <- get_all_p_vars(.GlobalEnv)
  if(print_output) cat("Found", length(p_vars), "p.* variables in global environment\n")
  if (length(p_vars) == 0) {
    # Always print errors
    cat("ERROR: No p.* variables found in global environment.\n")
    return(list())
  }
  
  # Generate random draws
  if(print_output) cat("Generating random draws...\n")
  random_draws <- generate_random_p_draws(print_output = print_output)
  
  if(print_output) cat("Generated", length(random_draws), "random draws\n")
  
  if (length(random_draws) > 0) {
    # Save original values
    original_values <- p_vars
    
    # Update global probability variables with random draws
    for (name in names(random_draws)) {
      assign(name, random_draws[[name]], envir = .GlobalEnv)
      if(print_output) cat("Assigned", name, "=", random_draws[[name]], "\n")
    }
    
    # Only restore original values if explicitly requested
    if (restore_originals) {
      for (name in names(original_values)) {
        assign(name, original_values[[name]], envir = .GlobalEnv)
      }
      if(print_output) cat("Restored original values\n")
    }
  } else {
    # Always print errors
    cat("ERROR: No random values generated\n")
  }
  
  # Return nothing by default
  if(print_output) return(random_draws) else return(invisible(NULL))
}

# Improved generate_gamma_params function to handle NaNs and warnings more robustly
generate_gamma_params <- function(val, print_output = FALSE) {
  # Ensure we have positive values to work with
  if(val <= 0) {
    if(print_output) cat("  Warning: Value", val, "is not positive, using fallback parameterization\n")
    return(list(shape=25, scale=val/20))  # Return reasonable parameters for small values
  }
  
  if (exists("gamma.parms.from.quantiles")) {
    # For cost values, we want relative quantiles around the value
    q1 <- val * 0.8  # Lower bound at 80% of value
    q2 <- val * 1.2  # Upper bound at 120% of value
    
    # First try: use suppressWarnings to prevent NaN warnings from bubbling up
    output <- tryCatch({
      suppressWarnings({
        params <- gamma.parms.from.quantiles(c(q1, q2), p=c(0.25, 0.75))
      })
      
      # Check if parameters are valid (not NA, not negative, not producing NaN)
      if(is.na(params$a) || is.na(params$b) || params$a <= 0 || params$b <= 0) {
        stop("Invalid gamma parameters detected")
      }
      
      # Try to use parameters in pgamma to see if they work
      test <- pgamma(q1, shape=params$a, rate=1/params$b)
      if(is.nan(test)) {
        stop("Parameters produce NaN in pgamma")
      }
      
      # Parameters are good, return them
      return(params)
    }, error = function(e) {
      # If that fails, use a simpler method with mean and CV
      if(print_output) cat("  Using direct parameterization for value:", val, "\n")
      # Use 0.2 as coefficient of variation (SD/mean)
      cv <- 0.2
      shape <- 1/(cv^2)  # Alpha parameter
      scale <- val/shape  # Beta parameter (scale)
      return(list(shape=shape, scale=scale))
    })
    
    return(output)
  } else {
    warning("gamma.parms.from.quantiles function not found")
    # Fallback to basic parameterization
    cv <- 0.2
    shape <- 1/(cv^2)
    scale <- val/shape
    return(list(shape=shape, scale=scale))
  }
}

# Fixed generate_random_c_draws with debug output and better vector handling
generate_random_c_draws <- function(print_output = FALSE) {
  c_vars <- get_all_c_vars(.GlobalEnv)
  if(print_output) cat("Processing", length(c_vars), "cost variables\n")
  random_draws <- list()
  
  for(name in names(c_vars)) {
    c_var <- c_vars[[name]]
    
    # Print variable information only if requested
    if(print_output) {
      if(is.vector(c_var) && length(c_var) > 1) {
        cat("Variable:", name, "is a vector with", length(c_var), "elements:", head(c_var), "...\n")
      } else {
        cat("Variable:", name, "=", c_var, "\n")
      }
    }
    
    if(is.vector(c_var) && length(c_var) > 1) {
      # Process vector costs element by element
      result_vector <- numeric(length(c_var))
      for(i in seq_along(c_var)) {
        val <- c_var[i]
        parms <- generate_gamma_params(val, print_output)
        if(is.na(parms)[1]) {
          if(print_output) cat("  Element", i, "value", val, "is invalid, keeping original\n")
          result_vector[i] <- val
        } else {
          a <- if(!is.null(parms$shape)) parms$shape else parms[1]
          scale <- if(!is.null(parms$scale)) parms$scale else parms[2]
          result_vector[i] <- rgamma(1, shape = a, scale = scale)
          if(print_output) cat("  Element", i, "value", val, "-> random draw:", result_vector[i], "\n")
        }
      }
      random_draws[[name]] <- result_vector
      if(print_output) cat("  Saved vector result for", name, "\n")
    } else {
      # Handle single value
      val <- c_var
      parms <- generate_gamma_params(val, print_output)
      if(is.na(parms)[1]) {
        if(print_output) cat("  Value", val, "is invalid, keeping original\n")
        random_draws[[name]] <- val
      } else {
        a <- if(!is.null(parms$shape)) parms$shape else parms[1]
        scale <- if(!is.null(parms$scale)) parms$scale else parms[2]
        draw <- rgamma(1, shape = a, scale = scale)
        random_draws[[name]] <- draw
        if(print_output) {
          cat("  Generated parameters for", name, ": a =", a, "scale =", scale, "\n")
          cat("  Random draw:", draw, "\n")
        }
      }
    }
  }
  
  if(print_output) cat("Final random_draws list contains", length(random_draws), "items\n")
  return(random_draws)
}

# Function to generate new random cost values and update global environment
generate_new_random_c_values <- function(return_originals = FALSE, print_output = FALSE) {
  # Get all cost variables from global environment
  c_vars <- get_all_c_vars(.GlobalEnv)
  if(print_output) {
    cat("Found", length(c_vars), "c.* variables in global environment\n")
    cat("c.* variables found:\n")
    print(names(c_vars))
  }
  
  if (length(c_vars) == 0) {
    # Always print errors
    cat("ERROR: No c.* variables found in global environment.\n")
    return(list())
  }
  
  # Generate random draws
  if(print_output) cat("Generating random c draws...\n")
  random_c_draws <- generate_random_c_draws(print_output = print_output)
  
  if(print_output) cat("Generated", length(random_c_draws), "random c draws\n")
  
  # Print info about generated values
  if (length(random_c_draws) > 0) {
    if(print_output) {
      cat("Random c draws:\n")
      for (name in names(random_c_draws)) {
        cat(name, "=", random_c_draws[[name]], "\n")
      }
    }
    
    # Update global environment
    for (name in names(random_c_draws)) {
      assign(name, random_c_draws[[name]], envir = .GlobalEnv)
    }
    
    # Return the random values if requested
    if(print_output) return(random_c_draws) else return(invisible(NULL))
  } else {
    # Always print errors
    cat("ERROR: No random c values generated\n")
    return(invisible(NULL))
  }
}

# Improved generate_random_u_draws function with full debugging support
generate_random_u_draws <- function(print_output = FALSE) {
  u_vars <- get_all_u_vars(.GlobalEnv)  # Use global environment like the other functions
  if(print_output) cat("Processing", length(u_vars), "utility variables\n")
  random_draws <- list()
  
  for(name in names(u_vars)) {
    u_val <- u_vars[[name]]
    
    # Print variable information only if requested
    if(print_output) {
      if(is.vector(u_val) && length(u_val) > 1) {
        cat("Variable:", name, "is a vector with", length(u_val), "elements:", head(u_val), "...\n")
      } else {
        cat("Variable:", name, "=", u_val, "\n")
      }
    }
    
    if(is.vector(u_val) && length(u_val) > 1) {
      # Process vector utilities element by element
      result_vector <- numeric(length(u_val))
      for(i in seq_along(u_val)) {
        val <- u_val[i]
        if(!is.numeric(val)) {
          if(print_output) cat("  Element", i, "value", val, "is not numeric, keeping original\n")
          result_vector[i] <- val
        } else if(val == 0) {
          if(print_output) cat("  Element", i, "value", val, "is zero, keeping as zero\n")
          result_vector[i] <- 0
        } else if(val <= 0 || val >= 1) {
          if(print_output) cat("  Element", i, "value", val, "is outside (0,1), keeping original\n")
          result_vector[i] <- val
        } else {
          alpha <- val * 10
          beta <- (1 - val) * 10
          result_vector[i] <- rbeta(1, alpha, beta)
          if(print_output) cat("  Element", i, "value", val, "-> random draw:", result_vector[i], "\n")
        }
      }
      random_draws[[name]] <- result_vector
      if(print_output) cat("  Saved vector result for", name, "\n")
    } else {
      # Handle single value
      if(!is.numeric(u_val)) {
        if(print_output) cat("  Value", u_val, "is not numeric, keeping original\n")
        random_draws[[name]] <- u_val
      } else if(u_val == 0) {
        if(print_output) cat("  Value", u_val, "is zero, keeping as zero\n")
        random_draws[[name]] <- 0
      } else if(u_val <= 0 || u_val >= 1) {
        if(print_output) cat("  Value", u_val, "is outside (0,1), keeping original\n")
        random_draws[[name]] <- u_val
      } else {
        alpha <- u_val * 10
        beta <- (1 - u_val) * 10
        draw <- rbeta(1, alpha, beta)
        random_draws[[name]] <- draw
        if(print_output) {
          cat("  Generated parameters: alpha =", alpha, "beta =", beta, "\n")
          cat("  Random draw:", draw, "\n")
        }
      }
    }
  }
  
  if(print_output) cat("Final random_draws list contains", length(random_draws), "items\n")
  return(random_draws)
}

# Function to generate new random utility values and update global environment
generate_new_random_u_values <- function(return_originals = FALSE, print_output = FALSE) {
  # Get all utility variables from global environment
  u_vars <- get_all_u_vars(.GlobalEnv)
  if(print_output) {
    cat("Found", length(u_vars), "u.* variables in global environment\n")
    cat("u.* variables found:\n")
    print(names(u_vars))
  }
  
  if (length(u_vars) == 0) {
    # Always print errors
    cat("ERROR: No u.* variables found in global environment.\n")
    return(list())
  }
  
  # Generate random draws
  if(print_output) cat("Generating random u draws...\n")
  random_u_draws <- generate_random_u_draws(print_output = print_output)
  
  if(print_output) cat("Generated", length(random_u_draws), "random u draws\n")
  
  # Print info about generated values
  if (length(random_u_draws) > 0) {
    if(print_output) {
      cat("Random u draws:\n")
      for (name in names(random_u_draws)) {
        cat(name, "=", random_u_draws[[name]], "\n")
      }
    }
    
    # Update global environment
    for (name in names(random_u_draws)) {
      assign(name, random_u_draws[[name]], envir = .GlobalEnv)
      if(print_output) cat("Assigned new", name, "=", random_u_draws[[name]], "\n")
    }
    
    # Return the random values if requested
    if(print_output) return(random_u_draws) else return(invisible(NULL))
  } else {
    # Always print errors
    cat("ERROR: No random u values generated\n")
    return(invisible(NULL))
  }
}

# Function to initialize model parameters in the global environment
initialize_model_parameters <- function(force = TRUE) {
  # Check if parameters already exist
  if (!force && length(ls(pattern = "^p\\.", envir = .GlobalEnv)) > 0 && 
      length(ls(pattern = "^c\\.", envir = .GlobalEnv)) > 0 && 
      length(ls(pattern = "^u\\.", envir = .GlobalEnv)) > 0) {
    cat("Model parameters already exist in global environment.\n")
    cat("To force re-initialization, call with force=TRUE\n")
    return(FALSE)
  }
  
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

  return(TRUE)
}

log_simulation_cycle_results <- function(sim_results, n.i, v.n, output_json) {
  # sim_results: list of simulation results, can be single run or multiple replications
  # n.i: number of individuals
  # v.n: vector of state names
  # output_json: file to append JSON lines to
  
  # Remove file if it exists to start fresh
  if (file.exists(output_json)) file.remove(output_json)
  
  # Create a container for all replications
  all_replications <- list()
  
  # Check if we have a list of replications or a single run
  is_multiple_replications <- !is.null(names(sim_results)) && "all_results" %in% names(sim_results)
  
  if (is_multiple_replications) {
    replications <- sim_results$all_results
    n_reps <- length(replications)
  } else {
    # Single run case - treat as one replication
    replications <- list(sim_results)
    n_reps <- 1
  }
  
  # Process each replication
  for (rep_idx in 1:n_reps) {
    current_rep <- replications[[rep_idx]]
    
    # Get matrices from current replication for both ESD and Surgery
    m.M_esd <- current_rep$esd$m.M
    m.C_esd <- current_rep$esd$m.C
    m.E_esd <- current_rep$esd$m.E
    
    m.M_surgery <- current_rep$surgery$m.M
    m.C_surgery <- current_rep$surgery$m.C
    m.E_surgery <- current_rep$surgery$m.E
    
    # Divide QALYs by 4 for annualization
    m.E_esd <- m.E_esd / 4
    m.E_surgery <- m.E_surgery / 4
    
    n.t <- ncol(m.M_esd) - 1
    rep_cycles <- list()
    
    # Process each cycle
    for (cycle in 1:(n.t + 1)) {
      cycle_data <- list()
      
      # First, process ESD data
      # State distribution for ESD
      esd_state_counts <- as.list(table(factor(m.M_esd[, cycle], levels = v.n)))
      names(esd_state_counts) <- v.n
      
      # Costs and QALYs up to this cycle for ESD
      esd_costs <- rowSums(m.C_esd[, 1:cycle, drop = FALSE])
      esd_qalys <- rowSums(m.E_esd[, 1:cycle, drop = FALSE])
      
      # Surveillance costs for ESD (sum costs for patients in S1 and S2 states)
      esd_surv_costs <- numeric(n.i)
      for (i in 1:n.i) {
        # Sum costs only for cycles where patient was in S1 or S2
        surv_cycles_esd <- which(m.M_esd[i, 1:cycle] %in% c("S1", "S2"))
        if (length(surv_cycles_esd) > 0) {
          esd_surv_costs[i] <- sum(m.C_esd[i, surv_cycles_esd])
        }
      }
      
      # Statistics for ESD
      esd_data <- list(
        state_counts = esd_state_counts,
        cost = list(
          total = list(avg = mean(esd_costs), std = sd(esd_costs)),
          surveillance = list(avg = mean(esd_surv_costs), std = sd(esd_surv_costs))
        ),
        qaly = list(avg = mean(esd_qalys), std = sd(esd_qalys))
      )
      
      # Second, process Surgery data
      # State distribution for Surgery
      surgery_state_counts <- as.list(table(factor(m.M_surgery[, cycle], levels = v.n)))
      names(surgery_state_counts) <- v.n
      
      # Costs and QALYs up to this cycle for Surgery
      surgery_costs <- rowSums(m.C_surgery[, 1:cycle, drop = FALSE])
      surgery_qalys <- rowSums(m.E_surgery[, 1:cycle, drop = FALSE])
      
      # Surveillance costs for Surgery (sum costs for patients in S1 and S2 states)
      surgery_surv_costs <- numeric(n.i)
      for (i in 1:n.i) {
        # Sum costs only for cycles where patient was in S1 or S2
        surv_cycles_surgery <- which(m.M_surgery[i, 1:cycle] %in% c("S1", "S2"))
        if (length(surv_cycles_surgery) > 0) {
          surgery_surv_costs[i] <- sum(m.C_surgery[i, surv_cycles_surgery])
        }
      }
      
      # Statistics for Surgery
      surgery_data <- list(
        state_counts = surgery_state_counts,
        cost = list(
          total = list(avg = mean(surgery_costs), std = sd(surgery_costs)),
          surveillance = list(avg = mean(surgery_surv_costs), std = sd(surgery_surv_costs))
        ),
        qaly = list(avg = mean(surgery_qalys), std = sd(surgery_qalys))
      )
      
      # Calculate ICER at this cycle
      delta_cost <- mean(surgery_costs) - mean(esd_costs)
      delta_qaly <- mean(surgery_qalys) - mean(esd_qalys)
      icer <- ifelse(delta_qaly != 0, delta_cost / delta_qaly, NA)
      
      # Add both ESD and Surgery data to cycle data
      cycle_data <- list(
        esd = esd_data,
        surgery = surgery_data,
        icer = list(avg = icer, std = NA)
      )
      
      # Add cycle data to rep_cycles
      rep_cycles[[as.character(cycle)]] <- cycle_data
    }
    
    # Add this replication to the results
    all_replications[[as.character(rep_idx)]] <- rep_cycles
  }
  
  # Write the entire result to one JSON file
  json_result <- jsonlite::toJSON(all_replications, auto_unbox = TRUE, pretty = TRUE)
  write(json_result, file = output_json)
}