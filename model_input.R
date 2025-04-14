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