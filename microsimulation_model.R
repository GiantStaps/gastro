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
#rm(list = ls())  # remove any variables in R's memory 

##################################### Model input #########################################
# Model input
n.i   <- 10000                 # number of simulated individuals
n.t   <- 20                    # time horizon, 20 cycles
v.n   <- c("H1","H2","S1",     # the model states: PreESD (H1), PreSurgery (H2), Surveillance-ESD (S1),
           "S2","P","D")       # Surveillance-Surgery(S2), Palliative (P), Death(D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H1", n.i)        # ESD model: everyone begins in PreESD
v.M_2 <- rep("H2", n.i)        # Surgery model: everyone begins in PreSurgery
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("ESD", "Surgery")   # store the strategy names

# Transition probabilities (per cycle)
p.H1D    <- 0.0050              # 3 month mortality rate after ESD
p.H1H2   <- 0.0886              # probability of switching to surgery after ESD
p.H1S1   <- 1 - p.H1D - p.H1H2  # probability of staying in surveillance after ESD

p.H2D    <- 0.0463             # short-term mortality rate after surgery
p.H2P    <- 0.0082             # probability of clinical failure after ESD
p.H2S2   <- 1 - p.H2P - p.H2D  # probability of staying in surveillance after surgery

# Pre-calculate S1 and S2 time-dependent transition probabilities
p.H1S1 <- 1 - p.H1D - p.H1H2  # probability of staying in surveillance after ESD
p.H2S2 <- 1 - p.H2P - p.H2D  # probability of staying in surveillance after surgery

# S1 state transition probabilities
p.S1H2 <- c(0,
            0.03,
            0.03,
            0.0249,
            0.0249,
            0.0249,
            0.0249,
            0.02)  # probability of local recurrence after ESD

p.S1P  <- rep(0, 8)     # probability of distant recurrence after ESD

p.S1D  <- c(0.002,
            0.0034,
            0.0034,
            0.0028,
            0.0028,
            0.0028,
            0.0028,
            0.0023)  # mortality rate for ESD follow-up

# Calculate p.S1S1 for each time point
p.S1S1 <- 1 - p.S1H2 - p.S1P - p.S1D  # probability of staying in S1

# S2 state transition probabilities
p.S2P  <- c(0.0000,
            0.0064,
            0.0064,
            0.0064,
            0.0038,
            0.0038,
            0.0038,
            0.0025)    # local + distant recurrence for surgery

p.S2D  <- c(0.0058,
            0.0058,
            0.0044,
            0.0044,
            0.0026,
            0.0026,
            0.0026,
            0.0018)    # mortality rate for surgery

# Calculate p.S2S2 for each time point
p.S2S2 <- 1 - p.S2P - p.S2D    # probability of staying in S2

# Palliative care probabilities
p.PD <- 0.5  # assume patient dies with 50% chance when in palliative care
p.PP <- 1 - p.PD  # probability of staying in palliative care


# Cost and utility inputs 
c.H1     <- 18347.00                 # cost of ESD
c.H2     <- 80000.00                # cost of surgery
c.S1     <- c(4137.50, 
2068.75, 
2068.75 ,
1443.75 ,
1443.75 ,
1443.75 ,
1443.75 ,
1443.75 
)                 # cost of surveillance-ESD 
c.S2     <- c(3275.00, 
 1637.50, 
 1637.50, 
 818.75, 
 818.75, 
 818.75, 
 818.75, 
 818.75)

c.P      <- 9600                  # cost of palliative

u.H1     <- 0.7                 # utility when PreESD 
u.H2     <- 0.7                 # utility when PreSurgery 
u.S1     <- c(0.73,
0.75,
0.75,
0.76,
0.76,
0.76,
0.76,
0.76
)                   # utility when surveillance-ESD
u.S2     <-  c(0.73,
0.75,
0.75,
0.76,
0.76,
0.76,
0.76,
0.76
)                   # utility when surveillance-surgery
u.P      <- 0.5                 # utility when palliative
u.D      <- 0                   # utility when dead

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


##################################### Run the simulation ##################################
# Run a smaller simulation for debugging/inspection purposes
n.i.debug <- 1000  # Use fewer individuals for debugging

# Run simulation with verbose output
cat("\n\nRunning ESD simulation:\n")
sim_esd_debug <- MicroSim(v.M_1[1:n.i.debug], n.i.debug, n.t, v.n, d.c, d.e, verbose = TRUE)

# Print ESD costs and utilities
cat("\n============== ESD Costs and Utilities ==============\n")
cat("Costs:\n")
cat("H1 (ESD procedure cost):", c.H1, "\n")
cat("S1 surveillance costs:", paste(c.S1, collapse=", "), "\n")
cat("Palliative care cost:", c.P, "\n")

cat("\nUtilities:\n")
cat("H1 (PreESD):", u.H1, "\n")
cat("S1 surveillance utilities:", paste(u.S1, collapse=", "), "\n")
cat("Palliative care utility:", u.P, "\n")
cat("Death utility:", u.D, "\n")

cat("\n============== Surgery Costs and Utilities ==============\n")
cat("Costs:\n")
cat("H2 (Surgery procedure cost):", c.H2, "\n")
cat("S2 surveillance costs:", paste(c.S2, collapse=", "), "\n")
cat("Palliative care cost:", c.P, "\n")

cat("\nUtilities:\n")
cat("H2 (PreSurgery):", u.H2, "\n")
cat("S2 surveillance utilities:", paste(u.S2, collapse=", "), "\n")
cat("Palliative care utility:", u.P, "\n")
cat("Death utility:", u.D, "\n")
cat("\n================================================\n\n")

cat("\n\nRunning Surgery simulation:\n")
sim_surgery_debug <- MicroSim(v.M_2[1:n.i.debug], n.i.debug, n.t, v.n, d.c, d.e, verbose = TRUE)

# Run the full simulations for cost-effectiveness analysis (no verbose output to avoid clutter)
sim_esd <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, verbose = FALSE)
sim_surgery <- MicroSim(v.M_2, n.i, n.t, v.n, d.c, d.e, verbose = FALSE)

################################# Cost-effectiveness analysis #############################

# store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v.C  <- c(sim_esd$tc_hat, sim_surgery $tc_hat) 
se.C <- c(sd(sim_esd$tc), sd(sim_surgery $tc)) / sqrt(n.i)
# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)
v.E  <- c(sim_esd$te_hat, sim_surgery $te_hat)
se.E <- c(sd(sim_esd$te), sd(sim_surgery $te)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # calculate incremental costs
delta.E <- v.E[2] - v.E[1]                   # calculate incremental QALYs
se.delta.E <- sd(sim_surgery $te - sim_esd$te) / sqrt(n.i) # Monte Carlo squared error (MCSE) of incremental costs
se.delta.C <- sd(sim_surgery $tc - sim_esd$tc) / sqrt(n.i) # Monte Carlo squared error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # calculate the ICER
results <- c(delta.C, delta.E, ICER)         # store the values in a new variable

# Create full incremental cost-effectiveness analysis table
table_micro <- data.frame(
  c(round(v.C, 0),  ""),           # costs per arm
  c(round(se.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # health outcomes per arm
  c(round(se.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # incremental costs
  c("", round(se.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # incremental QALYs 
  c("", round(se.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
rownames(table_micro) <- c(v.Trt, "* are MCSE values")  # name the rows
colnames(table_micro) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro  # print the table
print(table_micro)

# NEW: Print out QALYs and ICER for both procedures
cat("\n============== QALYs and ICER ==============\n")
cat("ESD QALYs:", round(sim_esd$te_hat, 3), "\n")
cat("Surgery QALYs:", round(sim_surgery$te_hat, 3), "\n")
cat("ICER:", round(ICER, 0), "\n")