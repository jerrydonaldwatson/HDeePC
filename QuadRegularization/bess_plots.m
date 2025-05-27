%% BESS Simulation (Example 1, config)
% Jeremy Watson, 2025
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
%
% Run this file to generate the comparison figure
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.
clear
close all
clc

load deepc_sq_reg

u_sim_dpc = u_sim;
y_sim_dpc = true_y_sim;

load hdeepc_sq_reg

u_sim_hdpc = u_sim;
y_sim_hdpc = true_y_sim;

figure(2)
subplot(3,1,1)
plot(1:k,u_sim_dpc(1,1:k), 1:k,u_sim_hdpc(1,1:k))
title('BESS current')
ylabel("Current (A)")
legend('DeePC', 'HDeePC', 'Location', 'best')

subplot(3,1,2)
plot(1:k,y_sim_dpc(1,1:k), 1:k,y_sim_hdpc(1,1:k))
title('Node voltage deviation')
ylabel("Voltage (V)")
legend('DeePC', 'HDeePC', 'Location', 'best')

subplot(3,1,3)
plot(1:k,y_sim_dpc(2,1:k), 1:k,y_sim_hdpc(2,1:k), 'LineWidth', 2)
title('BESS State-of-Charge (SoC)')
xlabel("Time step")
legend('DeePC', 'HDeePC', 'Location', 'best')