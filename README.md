# Microsimulation: ESD vs. Surgery Procedure and Recovery

This project uses **R version 4.5** but should be compatible with earlier R versions. The simulation is based on microsimulation methods from [DARTH-git's Microsimulation Tutorial](https://github.com/DARTH-git/Microsimulation-tutorial).

## Objective

The simulation models and compares patient experiences undergoing Endoscopic Submucosal Dissection (ESD) versus traditional surgical procedures, focusing on procedural outcomes and recovery trajectories.

## Features

- Simulates cohorts of patients simultaneously
- Compares health outcomes and cost-effectiveness of ESD and surgery
- Supports sensitivity analysis to evaluate parameter uncertainty
- Includes a Shiny GUI for easy simulation execution and visualization

## Running the Simulation

To run the simulation:

```r
# Make sure you're already in the gastro directory:
shiny::runApp("run_simulation.R")
```

The `run_simulation.R` script will automatically launch the Shiny GUI, allowing you to:
1. Configure simulation parameters (number of replications, individuals, etc.)
2. Run baseline or sensitivity analysis simulations
3. View real-time progress and results
4. Report state distribution for each cycle in results/ as a json file with date+time as its name
5. Explore visualizations of the simulation outcomes

## Dependencies

- R (4.5 recommended, but compatible with earlier versions)
- shiny
- dplyr
- ggplot2
- tidyr
- jsonlite
- future
- promises
- shinyjs

## Reference

- [Microsimulation Tutorial by DARTH-git](https://github.com/DARTH-git/Microsimulation-tutorial)
