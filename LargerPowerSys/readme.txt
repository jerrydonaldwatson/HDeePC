This folder contains the files needed to run Case Study 3 (Larger power system).
Jeremy Watson, 2025

This code accompanies our paper: "Hybrid Data-Enabled Predictive
Control: Incorporating model knowledge into the DeePC"
Example 3 (Parameter variation resulting in model mismatch)

Contact: jeremy.watson@canterbury.ac.nz in case of any queries.

Acknowledgement:
power_grid_uXX.mat and the example setting come from
U. Wasekar, J. Watson "Monte-Carlo analysis of interlinking converter
modelling and control in hybrid AC/DC networks", PMAPS2024.

Notes:
XX = level of parameter variation uncertainty
power_grid_u00.mat is therefore the base model
Files for both multiplicative and low-rank additive perturbation are included as shown in Table 3.

Please note, Mosek is required to run these cases and version 9.19 was used via an academic license. 
The code has been tested in MATLAB 2021b and 2023b. Other versions of MATLAB are likely to be able to run this, but not guaranteed.
In the case of any issues, please contact the author above. 

