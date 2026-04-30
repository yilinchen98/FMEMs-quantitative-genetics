# FMEMs-quantitative-genetics

This repository stores the R code for two models in **A functional mixed-effects model for the joint analysis of amplitude and phase variations in function-valued traits** by Yilin Chen, Davide Pigoli, John A.D. Aston and Patrick A. Carter. 

## Repository structure

- `PrepareR.R`: contains packages and functions required to fit the proposed genetic mixed-effect model.
- `AlignedFMEM/`
1. `Section4_simulation_study/`: contains R code to simulate datasets used in the simulation study.

    - `Sec4_1_Simulation_design.R`: This code simulate aligned curve data used in Section 4.1. 
    - `Sec4_2_Simulation_results.R`: This code reproduces results for fitting the simulated data with the proposed genetic mixed-effects model displayed in Figure 1 - Figure 3.
    - `Sec4_3_Bootstrap_simulation_results.R`: This code reproduces results for the bootstrap simulation and simultaneous confidence bands disaplayed in Figure 4 and Figure 5.
    - `Simulation_boostrapping_process.R`:This code give an example of fit aligned FMEM to one of the simulated dataset and then carry on bootstrapping.
    - `sim_results/`: contains simulated datasets, estimated fixed effect and covaraince functions and boostrap samples of the simulations.

2.
