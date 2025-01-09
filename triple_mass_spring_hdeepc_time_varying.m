%% Triple-Mass Simulation (Example 2)
% Jeremy Watson, 2025
%
% Acknowledgement:
% sys_dc.mat and the example setting come from
% https://github.com/4flixt/DeePC_Perspective/tree/main
% F. Fiedler "On the relationship between data-enabled predictive control
% and subspace predictive control", ECC 2021.
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
% This is the time-varying version
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
% Can be adjusted to Gaussian noise with \sigma = noiseM by changing rand()
% to randn() in appropriate places.

% Define System Matrices
load sys_dc.mat
A = A_dc;
B = B_dc;
C = C;
D = D;

%% DeePC Controller

% Adjust the number of unknown states n_u:
% n_u = 0 => MPC
% 0 < n_u < n => HDeePC
% n_u = n => DeePC
% Also ensure that the p_k is set correctly
m = size(B,2);          % Number of inputs
p = size(C,1);          % Number of outputs
p_k = 0;                % Number of outputs for which the equation is known
p_u = p - p_k;          % Number of outputs for which the equation is unknown
n = size(A,1);          % Number of states
n_u = 8;                % Number of unknown states
n_k = n - n_u;          % Number of known states
n_ti = 2;               % Number of time-invariant states (starting from x_1, x_2 etc.)

% Parameters
Tini = 4;       % Past data for learning the system
N = 20;         % Horizon
T = 150;        % Length of initial simulation to gather data

% Cost function weights
Q = 1; % All states and inputs are currently weighted the same
R = 0.1;
Qexpanded = diag(Q * ones(1,N*p_u));
Rexpanded = diag(R * ones(1,N*m));

% Regularization parameters
lg = 1e0;   % prevents overfitting
ly = 1e6; % helps deal with noise
lu = 1e6;

% Useful dimensions
L = Tini + N;
G = T-Tini-N+1;

% Targets to try to reach
r = zeros(N*p_u,1); % Used in optimisation problem
r_k = zeros(N*p_k,1); % Used in optimisation problem
r_c = zeros(p,1); % Used only to compute cost

% HDeePC partitioning
B_c = B(1:n_u,:);
B_k = B(n_u+1:end,:);
C_u = C(1:p_u,1:n_u);
C_f = C(1:p_u,n_u+1:end);
C_c = C(p_u+1:end,1:n_u);
C_k = C(p_u+1:end,n_u+1:end);
D_k = D(p_u+1:end,:);

C_y = C_c(1:p_k,1:p_u); % remember to adjust if y_u =/= x_u

% Adjust A with time k
for k = 1:max(lengthSim,T+1)
    A_tv(:,:,k) = A;
    A_tv(n_ti+1:end,:,k) = A(n_ti+1:end,:) + A(n_ti+1:end,:).*randn(size(A(n_ti+1:end,:))).*timeVarStd;

    A_u(:,:,k) = A_tv(1:n_u,1:n_u,k);
    A_f(:,:,k) = A_tv(1:n_u,n_u+1:end,k);
    A_c(:,:,k) = A_tv(n_u+1:end,1:n_u,k);
    A_k(:,:,k) = A_tv(n_u+1:end,n_u+1:end,k);

    A_y(:,:,k) = A_c(1:n_k,1:p_u,k); % remember to adjust if y_u =/= x_u

end

%% Initial simulation to gather data
u_d = 6 .* (rand(m, T+1) - 0.5);
y_d = zeros(p_u, T); % p_u outputs for data-based part of model
y_d_all = zeros(p, T); % used to record all outputs for info
x_d(:,1) = rand(n,1);
for k = 1:T+1
    x_d(:, k+1) = A_tv(:,:,k)*x_d(:, k) + B*u_d(:, k);

    noise = noiseM*2*(rand(p,1) - 0.5);
    y_d(:, k) = [C_u C_f]*x_d(:,k) + noise(1:p_u,1);
    y_d_all(:, k) = C*x_d(:,k) + noise;
end

% Reformat data matrices to build Hankel matrices
for i = 1:Tini*m+N*m
    for j = 1:G
        I = rem(i,m);
        if rem(i,m) == 0
            I = m;
        end
        H_u(i,j) = u_d(I, floor((i-1)/m + j));
    end
end

for i = 1:Tini*p_u+N*p_u
    for j = 1:G
        I = rem(i,p_u);
        if rem(i,p_u) == 0
            I = p_u;
        end
        H_y(i,j) = y_d(I, floor((i-1)/p_u + j));
    end
end

% Check rank (needs to be L)
if rank(H_u) >= L
    fprintf('Persistency of excitation satisfied!\n')
else
    fprintf('Persistency of excitation not satisfied!\n')
end

% Split Hankel matrices into past and future
if n_u > 0
    Up = H_u(1:Tini*m,:);
    Uf = H_u(Tini*m+1 : end,:);
    Yp = H_y(1:Tini*p_u,:);
    Yf = H_y(Tini*p_u+1 : end,:);

    % Initial data for learning the system
    uini = [];
    yini = [];
    for i = 1:Tini
        uini = [uini; u_d(:,size(u_d,2)-Tini+i)];
        yini = [yini; y_d(:,size(y_d,2)-Tini+i)];
    end
end

% Initial states
y_sim(1:p,1) = y_d_all(1:p,T + 1);
x(:, 1) = x_d(:,T + 2);
true_y_sim = y_sim; % true system dynamics without measurement noise, for plots

cvx_solver sdpt3

%% Main loop
cost_u(lengthSim) = 0;
cost_y(lengthSim) = 0;
tic
for k = 1:lengthSim

    % Solve optimisation problem using cvx
    cvx_begin quiet

    variable g(G,1)
    variable u(N*m,1);
    variable y_u(N*p_u,1);
    variable y_k(N*p_k,1);
    variable x_k(n_k,N+1);
    variable sigma_y(Tini*p_u,1);
    variable sigma_u(Tini*m,1);

    minimize quad_form(y_u-r,Qexpanded) + quad_form(y_k-r_k,Q*eye(N*p_k)) +  quad_form(u,Rexpanded) + quad_form(g,lg*eye(G)) + lu*norm(sigma_u,1) + ly*norm(sigma_y,1)%quad_form(sigma_u,lu*eye(Tini*m)) + quad_form(sigma_y,ly*eye(Tini*p_u)) % Complete cost function

    % Data-based part
    if n_u > 0
        [Up; Yp; Uf; Yf] * g == [(uini + sigma_u); (yini + sigma_y); u; y_u] % trajectories have to be consistent with initial data
    end

    % Constraints on u and (possibly) y
    for i = 1:m*N
        u(i) <= 0.7
        u(i) >= -0.7
    end

    % MPC part
    x_k(:,1) == x(n_u+1:end, k);
    for i = 1:N
        x_k(:,i+1) == A_y(:,:,k+i-1) * y_u((1:p_u)+(i-1)*p_u) + A_k(:,:,k+i-1) * x_k(:,i) + B_k * u((1:m)+(i-1)*m); % note that y_u(1:2) = x_u(1:2)
        y_k((1:p_k)+(i-1)*p_k) == C_y * y_u((1:p_u)+(i-1)*p_u) + C_k * x_k(:,i) + D_k * u((1:m)+(i-1)*m);
    end

    cvx_end

    % Assign outputs for next time-step (at present, recalculating every step...)
    if ~isnan(u(1))
        u_sim(:,k) = u(1:m,1);
    else
        fprintf('Optimization problem failed to solve \n')
        u_sim(:,k) = 0;
    end

    % Update system
    x(:, k+1) = A_tv(:,:,k) * x(:, k) + B * u_sim(:, k);
    y_sim(:, k) = C * x(:, k) + D * u_sim(:, k) + noiseM*2*(rand(p,1) - 0.5);
    true_y_sim(:, k) = C * x(:, k) + D * u_sim(:, k);

    % Update data vectors - new data at the end of the vector
    if n_u > 0
        uini(1:m,:) = [];
        yini(1:p_u,:) = [];
        uini = [uini; u_sim(:, k)];
        yini = [yini; y_sim(1:p_u, k)];
    end

    % Record costs
    cost_y(k) = Q * x(1:p, k+1)'*x(1:p, k+1); % Adjust for a different cost function
    cost_u(k) = R * u_sim(:, k)'*u_sim(:, k);
end
toc

% Total cost:
fprintf('Cost is %6.4f\n',sum(cost_y + cost_u));

% Plot DeePC input, predicted output, and reference
figure(1)
title("DeePC - triple-mass system")
subplot(3,1,1)
plot(1:k,u_sim(:,1:k))
legend('u_1', 'u_2')

subplot(3,1,2)
plot(1:k,true_y_sim(1:3,1:k))
legend('y_1', 'y_2', 'y_3')

subplot(3,1,3)
plot(1:k,cost_y + cost_u)
legend('Cost')