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
p.H1D    <- 0.005              # probability to die after ESD
p.H1H2   <- 0.05
p.H1S1   <- 0.945
p.H2D    <- 0.005              # probability to die after surgery
p.H2S2   <- 0.945         	      
p.H2P    <- 0.05
p.S1H2   <- 0.05
p.S1S1   <- 0.945 
p.S1D    <- 0.005 
p.S2S2   <- 0.945 
p.S2P    <- 0.05
p.S2D    <- 0.005
p.PP     <- 0.5
p.PD     <- 0.5
   

# Cost and utility inputs 
c.H1     <- 4000                # cost of ESD
c.H2     <- 4000                # cost of surgery
c.S1     <- 100                 # cost of surveillance-ESD 
c.S2     <- 100                 # cost of surveillance-surgery
c.P      <- 50                  # cost of palliative

u.H1     <- 0.8                 # utility when PreESD 
u.H2     <- 0.8                 # utility when PreSurgery 
u.S1     <- 1                   # utility when surveillance-ESD
u.S2     <- 1                   # utility when surveillance-surgery
u.P      <- 0.5                 # utility when palliative

##################################### Functions ###########################################

# The MicroSim function for the simple microsimulation of the 'Sick-Sicker' model keeps track of what happens to each individual during each cycle. 

MicroSim <- function(v.M, n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, seed = 1) {
  # Arguments:  
  # v.M:     vector of initial states for individuals 
  # n.i:     number of individuals
  # n.t:     total number of cycles to run the model
  # v.n:     vector of health state names
  # d.c:     discount rate for costs
  # d.e:     discount rate for health outcome (QALYs)
  # TR.out:  should the output include a microsimulation trace? (default is TRUE)
  # TS.out:  should the output include a matrix of transitions between states? (default is TRUE)
  # seed:    starting seed number for random number generator (default is 1)
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
  
  m.M[, 1] <- v.M                     # indicate the initial health state   
  
  for (i in 1:n.i) {
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m.C[i, 1] <- Costs(m.M[i, 1])  # estimate costs per individual for the initial health state 
    m.E[i, 1] <- Effs (m.M[i, 1])  # estimate QALYs per individual for the initial health state
    
    for (t in 1:n.t) {
      v.p <- Probs(m.M[i, t])           # calculate the transition probabilities at cycle t 
      
      m.M[i, t + 1] <- sample(v.n, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      m.C[i, t + 1] <- Costs(m.M[i, t + 1])   # estimate costs per individual during cycle t + 1
      m.E[i, t + 1] <- Effs( m.M[i, t + 1])   # estimate QALYs per individual during cycle t + 1
      
    } # close the loop for the time points 
    if(i/100 == round(i/100,0)) {          # display the progress of the simulation
      cat('\r', paste(i/n.i * 100, "% done", sep = " "))
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
  
  results <- list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  v.p.it <- rep(NA, n.s)     # create vector of state transition probabilities
  names(v.p.it) <- v.n       # name the vector
  
  # update v.p.it with the appropriate probabilities   
  v.p.it[M_it == "H1"]  <- c(0, p.H1H2, p.H1S1, 0, 0, p.H1D)                 # transition probabilities when H1
  v.p.it[M_it == "H2"]  <- c(0, 0, 0, p.H2S2, p.H2P, p.H2D)                  # transition probabilities when H2
  v.p.it[M_it == "S1"]  <- c(0, p.S1H2, p.S1S1, 0, 0, p.S1D)                 # transition probabilities when S1
  v.p.it[M_it == "S2"] <- c(0, 0, 0, p.S2S2, p.S2P, p.S2D)                   # transition probabilities when S2
  v.p.it[M_it == "P"]  <- c(0, 0, 0, 0, p.PP, p.PD)                          # transition probabilities when P
  v.p.it[M_it == "D"]  <- c(0, 0, 0, 0, 0, 1)                                # transition probabilities when dead   
  ifelse(sum(v.p.it) == 1, return(v.p.it), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}       


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function (M_it) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  c.it <- 0                                  # by default the cost for everyone is zero 
  c.it[M_it == "H1"]  <- c.H1                # update the cost if H1
  c.it[M_it == "H2"]  <- c.H2                # update the cost if H2
  c.it[M_it == "S1"]  <- c.S1                # update the cost if S2
  c.it[M_it == "S2"]  <- c.S2                # update the cost if S2
  c.it[M_it == "P"]   <- c.P                 # update the cost if P
  return(c.it)        		                   # return the costs
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, cl = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[M_it == "H1"]  <- u.H1    # update the utility if H1
  u.it[M_it == "H2"]  <- u.H2    # update the utility if H2
  u.it[M_it == "S1"] <- u.S1     # update the utility if S1
  u.it[M_it == "S2"] <- u.S2     # update the utility if S2
  u.it[M_it == "P"]  <- u.P      # update the utility if P
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}


##################################### Run the simulation ##################################
sim_esd  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e) # run for ESD
sim_surgery     <- MicroSim(v.M_2, n.i, n.t, v.n, d.c, d.e)  # run for surgery
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