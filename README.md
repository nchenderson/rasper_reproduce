# rasper_reproduce

Code for reproducing results in the RASPER paper.

# Simulation Code

## Simulation1

This folder contains two scripts for running the 
simulations described in the manuscript.

Simulation1aRun.R - code for running and saving simulation
                    results for Simulation 1(a)
Simulation1bRun.R - code for running and saving simulation
                    results for Simulation 1(b)

## Simulation2

This folder contains two scripts for running the 
simulations described in the manuscript.

Simulation2aRun.R - code for running and saving simulation
                    results for Simulation 2(a)
Simulation2bRun.R - code for running and saving simulation
                    results for Simulation 2(b)

# GenerateTablesFigures

This folder contains scripts needed to directly generate
the tables and giures presented in the manuscript.

# Data Analysis

* FitRiskModels.R - Code for estimating the ICI prostate cancer using 
                    various approaches: OLS, ridge regression, direct
                    transfer learning, angle transfer learning, 
                    and two RASPER methods.

# Additional R scripts

* ExternalFunctions.R - This script contains code to perform estimation
                        and leave-one-out cross-validation from non-RASPER 
                        methods including ridge regression, ATL,
                        and DTL (informative ridge regression)
