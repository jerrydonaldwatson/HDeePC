%% BESS Simulation (Example 1, config)
% Jeremy Watson, 2025
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
%
% Run this file first to set-up parameters for both bess_hdeepc.m and
% bess_deepc.m
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.
clear
close all
clc

rng(1)

%% System Set Up
% Simulation parameters
lengthSim = 1000;
ts = 1e-3;        % Simulation time-step
noiseM = 1e-6;    % Noise magnitude time-step
randnGenLoad = 1; % Random generation / load at bus magnitude std
        
% Constraints
dV = 20; % maximum voltage deviation (V)
Imax = 5; % maximum current from BESS (A)

% Define power system / BESS parameters
Cdc = 1e-3;
Rline = 2;
Rload = 50;
Lline = 5e-3; % Line inductance
tq = 1e-4; % 1 / SoC time constant

% Generate random generation/load
u_genload = randn(1,lengthSim+N)*randnGenLoad;

%% DeePC Controller
% Horizons etc
Tini = 50;      % Past data for learning the system
N = 10;         % Horizon
T = 200;        % Initial simulation to gather data

% Cost function weights
Q = diag([0.001, 5e4]);
R = 0.001;

save config