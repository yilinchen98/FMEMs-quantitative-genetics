# FMEMs-quantitative-genetics

This repository stores the R code for two models in **A functional mixed-effects model for the joint analysis of amplitude and phase variations in function-valued traits** by Yilin Chen, Davide Pigoli, John A.D. Aston and Patrick A. Carter. 

## Repository structure

- `TRFUN25PUP4.DAT`: the clean-up flour beetle growth curves dataset (miss data removed).
- `PrepareR.R`: contains packages and functions required to fit the proposed genetic mixed-effect model.
- `AlignedFMEM/`:
1. `Section4_simulation_study/`: contains R code to simulate datasets used in the simulation study.

    - `Sec4_1_Simulation_design.R`: This code simulates aligned curve data used in Section 4.1. 
    - `Sec4_2_Simulation_results.R`: This code reproduces results for fitting the simulated data with the proposed genetic mixed-effects model displayed in Figure 1 - Figure 3.
    - `Sec4_3_Bootstrap_simulation_results.R`: This code reproduces results for the bootstrap simulation and simultaneous confidence bands displayed in Figure 4 and Figure 5.
    - `Simulation_boostrapping_process.R`: This code provides an example of fitting aligned FMEM to one of the simulated datasets and the bootstrapping procedure.
    - `sim_results/`: contains simulated datasets, estimated fixed effect and covariance functions and bootstrap samples of the simulations.

2. `TC_alignedFMEM_application/`:

    - `Sec5_3_TC_alignedFMEM_results.R`: This code reproduces results for fitting the aligned flour beetle growth data with the proposed genetic mixed-effects model displayed in Figure 7 - Figure 11. Bootstrap samples stored in `TC_CG_fun_bs_samples.rds`

- `JointFMEM/`:

    - `Sec5_4_TC_JointFMEM_results.R`: This code reproduces results for fitting the flour beetle growth data with the proposed joint genetic mixed-effects model displayed in Figure 12 - Figure 21. 
    - `JointFMEM_results/`: contains the estimated covariance functions and the bootstrap samples from fitting the joint model with the flour beetle growth data.
