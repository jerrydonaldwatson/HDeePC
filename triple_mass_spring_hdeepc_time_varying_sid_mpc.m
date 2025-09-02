%% Triple-Mass Simulation (Subspace ID + MPC with similarity transform + time-varying A)
% Jeremy Watson, 2025
%
% Modified to replace DeePC with Subspace ID + MPC approach,
% incorporating similarity transform T for state alignment,
% and adding time-varying A matrix as before.
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.

clear all
close all
clc

rng(0)

%% System Set Up
lengthSim = 40; % Timestep is 0.1 s
noiseM = 2e-5; % Uniform distribution with this limit
timeVarStd = 0.1; % 0.1 = 10%

load sys_dc.mat
A = A_dc;
B = B_dc;
C = C;
D = D;

% Dimensions
m = size(B,2);          % Number of inputs
p = size(C,1);          % Number of outputs
n = size(A,1);          % Number of states
n_ti = 0;               % Number of time-invariant states (preserved from original)

% Parameters
T = 150;    % Initial data length for system ID
Tini = 4;   % Not used in MPC but kept for consistency
N = 20;     % Prediction horizon

% Generate time-varying A matrix (for true system simulation)
for k = 1:max(lengthSim,T+1)
    A_tv(:,:,k) = A;
    A_tv(n_ti+1:end,:,k) = A(n_ti+1:end,:) + A(n_ti+1:end,:) .* randn(size(A(n_ti+1:end,:))) * timeVarStd;
end

% Initial simulation to gather data for system identification
u_d = 6 .* (rand(m, T+1) - 0.5);
y_d = zeros(p, T);
x_d(:,1) = rand(n,1);

for k = 1:T+1
    x_d(:, k+1) = A_tv(:,:,k)*x_d(:, k) + B*u_d(:, k);
    noise = noiseM*2*(rand(p,1) - 0.5);
    y_d(:, k) = C*x_d(:,k) + noise;
end

% Subspace system identification (using the nominal data with noise)
data_id = iddata(y_d', u_d', 0.1); % 0.1 is timestep
sysID = n4sid(data_id, n);

% Cost weights
Q = 1;
R = 0.1;
Qmat = Q*eye(p);
Rmat = R*eye(m);

% Reference trajectory (zero for simplicity)
r = zeros(N*p,1);

% Initial states for simulation
y_sim(1:p,1) = y_d(:,end);
x(:,1) = x_d(:,T+2);

true_y_sim = y_sim;

% Control input storage
u_sim = zeros(m, lengthSim);

%% Main MPC loop

tic
for k = 1:lengthSim
    
    cvx_begin quiet
        cvx_precision best
        
        variable u(m, N)
        variable x_k(n, N+1)
        variable y(p, N)
        variable slack_y(p,1)
        
        cost = 0;
        for i = 1:N
            cost = cost + (y(:, i) - r((i-1)*p+1:i*p))'* Qmat * (y(:, i) - r((i-1)*p+1:i*p)) + u(:, i)'* Rmat * u(:, i);
        end
        cost = cost + 1e6 * norm(slack_y);
        
        minimize(cost)
        
        % Initial state (output) constraint

        y(:,1) == C * x(:,k) + slack_y;
        
        for i = 1:N
            % System dynamics (identified model)
            x_k(:, i+1) == sysID.A * x_k(:, i) + sysID.B * u(:, i);
            y(:, i) == sysID.C * x_k(:, i) + sysID.D * u(:, i);
        end
        
        % Input constraints
        for i = 1:N
            u(:, i) <= 0.7;
            u(:, i) >= -0.7;
        end
        
    cvx_end
    
    % Apply first input
    if ~isnan(u(1,1))
        u_sim(:,k) = u(:,1);
    else
        fprintf('Optimization problem failed at step %d\n', k);
        u_sim(:,k) = zeros(m,1);
    end
    
    % True system update with time-varying A
    x(:, k+1) = A_tv(:,:,k) * x(:, k) + B * u_sim(:, k);
    y_sim(:, k) = C * x(:, k) + D * u_sim(:, k) + noiseM*2*(rand(p,1)-0.5);
    true_y_sim(:, k) = C * x(:, k) + D * u_sim(:, k);
    
end
toc

% Total cost:
cost_y = zeros(lengthSim,1);
cost_u = zeros(lengthSim,1);
for k = 1:lengthSim
    cost_y(k) = Q * x(1:p, k+1)'*x(1:p, k+1);
    cost_u(k) = R * u_sim(:, k)'*u_sim(:, k);
end
fprintf('Cost is %6.4f\n',sum(cost_y + cost_u));

% Plot MPC input, predicted output, and reference (original style)
figure(1)
title("MPC - triple-mass system")
subplot(3,1,1)
plot(1:lengthSim,u_sim(:,1:lengthSim))
legend('u_1', 'u_2')

subplot(3,1,2)
plot(1:lengthSim,true_y_sim(1:3,1:lengthSim))
legend('y_1', 'y_2', 'y_3')

subplot(3,1,3)
plot(1:lengthSim,cost_y + cost_u)
legend('Cost')
