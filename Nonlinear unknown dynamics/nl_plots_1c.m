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

load deepc

u_sim_dpc = u_sim;
y_sim_dpc = true_y_sim;

load hdeepc

u_sim_hdpc = u_sim;
y_sim_hdpc = true_y_sim;

load bess_nl_sid_mpc.mat

u_sim_mpc = u_sim;
y_sim_mpc = true_y_sim;

figure(2)


subplot(4,1,1)
plot(1:k,u_sim_dpc(1,1:k), 1:k,u_sim_hdpc(1,1:k), 1:k,u_sim_mpc(1,1:k), 'LineWidth', 2)
title('u_1')
%legend('DeePC', 'HDeePC', 'Sys ID + MPC', 'Location', 'best')
set(gca,'FontSize',12) 
subplot(4,1,2)
plot(1:k,u_sim_dpc(2,1:k), 1:k,u_sim_hdpc(2,1:k), 1:k,u_sim_mpc(2,1:k), 'LineWidth', 2)
title('u_1')
%legend('DeePC', 'HDeePC', 'Sys ID + MPC', 'Location', 'best')
set(gca,'FontSize',12) 
subplot(4,1,3)
plot(1:k,y_sim_dpc(1,1:k), 1:k,y_sim_hdpc(1,1:k), 1:k,y_sim_mpc(1,1:k), 'LineWidth', 2)
title('x_1')
%legend('DeePC', 'HDeePC', 'Location', 'best')
set(gca,'FontSize',12) 
subplot(4,1,4)
plot(1:k,y_sim_dpc(2,1:k), 1:k,y_sim_hdpc(2,1:k), 1:k,y_sim_mpc(2,1:k), 'LineWidth', 2)
title('x_3')
legend('DeePC', 'HDeePC', 'Sys ID + MPC', 'Target', 'Location', 'best')
xlabel('Time-step')
set(gca,'FontSize',12) 

figure(3)


subplot(2,1,1)
plot(1:k,u_sim_dpc(:,1:k), 1:k,u_sim_hdpc(:,1:k), 1:k,u_sim_mpc(:,1:k), 'LineWidth', 2)
title('Inputs')
legend('u_1: DeePC', 'u_2: DeePC', 'u_1: HDeePC', 'u_2: HDeePC', 'u_1: ID + MPC', 'u_2: ID + MPC', 'Location', 'best')
set(gca,'FontSize',12) 
subplot(2,1,2)
plot(1:k,y_sim_dpc(:,1:k), 1:k,y_sim_hdpc(:,1:k), 1:k,y_sim_mpc(:,1:k), 'LineWidth', 2)
title('Outputs')
legend('y_1: DeePC', 'y_2: DeePC', 'y_1: HDeePC', 'y_2: HDeePC', 'y_1: ID + MPC', 'y_2: ID + MPC', 'Location', 'best')
set(gca,'FontSize',12) 
xlabel("Time-step")
