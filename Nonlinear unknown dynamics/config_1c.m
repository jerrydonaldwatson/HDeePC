%% Example 1c, config
% Jeremy Watson, 2025
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
%
% Run this file first to set-up parameters for hdeepc_nonlinear_1c.m 
% and deepc_nonlinear_1c.m and nl_subspaceid_mpc_1c.m
% Example 1c
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.
clear
close all
clc

rng(1)

%% System Set Up
% Simulation parameters
lengthSim = 100;
ts = 1e-3;        % Simulation time-step
noiseM = 1e-6;    % Noise magnitude time-step
        
% Constraints
dX1 = 20;
dX3 = 20;
Umax = 5; 

% Define power system / BESS parameters
tq = 1;  % 1 for Fig. 4
e = 0.9; % 0.9 for Fig. 4
a = 0.7; % 0.7 for Fig. 4
b = 0.1; % 0.2 for Fig. 4

B1 = [1 1;
    0 1;
    -ts*tq/e 0];

B2 = [1 1;
    0 1;
    -ts*tq*e 0];

C = [1, 0, 0;
    0 0 1];

D = [0 0;
    0 0];

%% (H)DeePC Controller
% Horizons etc
Tini = 10;      % Past data for learning the system
N = 20;         % Horizon
T = 200;        % Initial simulation to gather data
maxIter = 10;   % Maximum iterations to solve non-convex problem iteratively with cvx

% Regularization parameters
lg = 1e3;   % prevents overfitting
ly = 1e6; % helps deal with noise

%% Sys ID + MPC
useSlack = true;
lx = 1e4;
ly = 1e4;

% Cost function weights
Q = diag([1e3, 1e3]);
R = 1;

save config