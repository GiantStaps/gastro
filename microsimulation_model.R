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


##################################### Generate Random Input ##################################
# Source all R files in the current directory (except this file itself)
r_files <- list.files(path = "c:/Users/rma86/Desktop/gastro", pattern = "\\.r$", ignore.case = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("microsimulation_model\\.r$", r_files, ignore.case = TRUE)]
for (f in r_files) source(f)


# Get all c.* and u.* variables (for potential future use)
p_vars <- get_all_p_vars()
c_vars <- get_all_c_vars()
u_vars <- get_all_u_vars()

# # Print first few random probabilities to verify
# cat("Sample of random probability draws:\n")
# if (length(random_p_draws) > 0) {
#   sample_names <- head(names(random_p_draws), 5)
#   for (name in sample_names) {
#     cat(name, "=", round(random_p_draws[[name]], 4), "\n")
#   }
# }

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