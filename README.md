# HDeePC
Accompanying repository for "Hybrid Data-enabled Predictive Control: Incorporating model knowledge into the DeePC".

Two case studies are presented, a BESS in a DC microgrid (Case study 1) and a triple-mass system (Case study 2).
The former is original, while the latter is derived from https://github.com/4flixt/DeePC_Perspective. 
The first study (part a) shows the application of HDeePC to improve the numerical stability of a system with slow/fast time-scales,
where the BESS state-of-charge dynamics are known but the microgrid dynamics are not. 
Part b (in the Nonlinear HDeePC) folder shows the application of HDeePC to a simple nonlinear BESS system. 

The second study shows the application of HDeePC to improve computational performance and reduce cost for: 
a) varying numbers of known states and outputs, and b) a linear time-varying system. 
The code is designed such that adjusting the number of known states to zero yields DeePC, while adjusting the number of known states to n (the order of the system) yields MPC.

Jeremy Watson, 2025: jeremy.watson (at) canterbury.ac.nz

# Files:
bess_config.m -  Run this first to set the configuration for the BESS example (Case study 1).

bess_deepc.m  -  This is the DeePC simulation

bess_hdeepc.m -  This is the HDeePC simulation

bess_plots.m  -  Run this after running all the files above to see the Case study 1 figure. 

sys_dc.mat           - This is the system data (A,B,C,D) file for Case study 2. Reproduced from https://github.com/4flixt/DeePC_Perspective. 

triple_mass_hdeepc.m - This runs Case study 2 for varying numbers of known state / output equations. Adjust p_k and n_u to replicate Table I. 

triple_mass_hdeepc.m - This runs Case study 2 for a time-varying system. Adjust p_k and n_u to replicate Table II. 

Condensed Optimization Code - This implements Case study 1a using a condensed optimization formulation, where the only decision variables are g and the slack variable (sigma_y).

Nonlinear HDeePC - This implements Case study 1b.

QuadRegularization - This changes the regularization to a quadratic form and again compares DeePC and HDeePC. 
