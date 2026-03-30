# Data-Driven Real-Time Optimization for Compressor Networks

This repository contains the MATLAB/Simulink implementation of a data-driven Real-Time Optimization (RTO) framework for parallel and serial centrifugal compressor networks. The algorithm utilizes Gaussian Processes (GPs) and Genetic Algorithms (GAs) to minimize power consumption while strictly satisfying dynamic surge limits and demand tracking constraints.


To ensure maximum compatibility with Simulink, this repository uses a flat file structure. All scripts, models, and dependencies are located in the root directory.

* **Main S-Functions:** `parallel_cei.m`, `parallel_ucb.m`, `series_cei.m`, `series_ucb.m`
* **Acquisition Functions:** `objectiveFun_*.m` files containing the rigorous mathematical penalty formulations.
* **Optimizers:** `solveOptimization_*.m` files containing the Genetic Algorithm solvers.
* **Surrogate Models:** `make_gp_*.m` and `gp_predict.m` handling the Bayesian learning and predictions.
* **Simulink Models:** `*.slx` files representing the physical compressor networks and steady-state detectors.


## Quick Start
1. Open MATLAB and navigate to the root directory of this repository. (No need to configure MATLAB paths; everything is in one folder).
2. Open the desired Simulink `.slx` model.
3. Ensure the S-Function block in the model is pointing to the correct main script (e.g., `series_cei`).
4. Run the simulation. 
5. **Tuning:** Hyperparameters (like safety margins, tracking weights, and target tolerances) can be modified directly at the bottom of the main S-function scripts inside the `initParams` function.

## Performance Note: Anti-Surge Control (ASC)
The detailed Anti-Surge Control (ASC) blocks have been reintegrated into the Simulink models to provide high-fidelity dynamic responses near the surge limit. However, this significantly increases the computational load and overall simulation time. 

**For faster simulations (Rapid Prototyping):** You can easily bypass the full ASC dynamics. Simply delete or disconnect the ASC block in the Simulink model and replace it with the **Step block** (labeled for recycle valve opening) located immediately next to it. This provides a simplified, much faster approximation of the valve dynamics, which is highly recommended when tuning the optimization hyperparameters.
