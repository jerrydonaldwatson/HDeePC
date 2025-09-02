This folder contains the files needed to run Case Study 1c (nonlinear system with unknown nonlinear part).
Jeremy Watson, 2025

This code accompanies our paper: "Hybrid Data-Enabled Predictive
Control: Incorporating model knowledge into the DeePC"
Example 1c (Nonlinear system with unknown nonlinear dynamics)

Contact: jeremy.watson@canterbury.ac.nz in case of any queries.

Notes:
The files show DeePC, HDeePC and system ID  + MPC. 
The ****_config.m file should be run first to ensure identical configuration for all three methods.
After all three controllers have been run, ****_plots should be run to generate the plot (Fig. 4 in the paper). 
The system identification toolbox is required for n4sid (system ID). 

The code has been tested in MATLAB 2021b and 2023b. Other versions of MATLAB are likely to be able to run this, but this is not guaranteed.
In the case of any issues, please contact the author above. 
