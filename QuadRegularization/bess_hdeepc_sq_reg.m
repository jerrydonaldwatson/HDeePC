%% BESS Simulation (Example 1, HDeePC - using a quadratic regularisation term)
% Jeremy Watson, 2025
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
% This example is referred to in the discussion about the regularization of
% Exmaple 1. 
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.


clear
close all
clc

rng(1)

%% System Set Up
load config

% Define System Matrices
A_d = [1-ts*(1/Rload/Cdc) ts/Cdc;%, 0; % voltage
    -ts*(1/Lline) 1-ts*Rline/Lline];
%    0, 1]; % SoC

B_d = [ts/Cdc ts/Cdc;
    0 0];
%    -ts*tq 0];

C_d = [1 0];
%    0 1];

D_d = [0 0];
%    0 0];

A = [1-ts*(1/Rload/Cdc), ts/Cdc 0; % voltage
    -ts*(1/Lline) 1-ts*Rline/Lline 0;
    0, 0, 1]; % SoC

B = [ts/Cdc ts/Cdc; 
    0 0;
    -ts*tq 0];

C = [1, 0, 0;
    0 0 1];

D = [0 0;
    0 0];

%% DeePC
m = 2;          % Number of inputs
p = 1;          % Number of outputs
n = 3;          % Number of states

% for I = 1:N
%     Qexpanded(2*I-1,2*I-1) = Q(1,1);
%     Qexpanded(2*I,2*I) = Q(2,2);
% end

Qexpanded = diag(Q(1,1) * ones(1,N));
Q2expanded = diag(Q(2,2) * ones(1,N));
Rexpanded = diag(R * ones(1,N*m));

% Regularization parameters
lg = 1;   % prevents overfitting
ly = 1e6; % helps deal with noise

% Initial simulation to gather data
u_d = 10 .* (rand(m, T+1) - 0.5);
y_d = zeros(p, T);
x_d(1,1) = rand(1,1);
x_d(2,1) = 0;
SOC0 = rand(1);
for k = 1:T+1
    x_d(:, k+1) = A_d*x_d(:, k) + B_d*u_d(:, k);
    y_d(:, k) = C_d*x_d(:,k) + noiseM*randn(p,1);
end

% Useful dimensions
L = Tini + N;
G = T-Tini-N+1;

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

for i = 1:Tini*p+N*p
    for j = 1:G
        I = rem(i,p);
        if rem(i,p) == 0
            I = p;
        end
        H_y(i,j) = y_d(I, floor((i-1)/p + j));
    end
end

% Check rank (needs to be L)
if rank(H_u) >= L
    fprintf('Persistency of excitation okay!\n')
else
    fprintf('Persistency of excitation not okay!\n')
end

% Split Hankel matrices into past and future
Up = H_u(1:Tini*m,:);
Uf = H_u(Tini*m+1 : end,:);
Yp = H_y(1:Tini*p,:);
Yf = H_y(Tini*p+1 : end,:);

% Initial data for learning the system
uini = [];
yini = [];
for i = 1:Tini
    uini = [uini; u_d(:,size(u_d,2)-Tini+i)];
    yini = [yini; y_d(:,size(y_d,2)-Tini+i)];
end

% Trajectory to try to follow
% for I = 1:N
%     r(2*I-1,1) = 0;
%     r(2*I,1) = 0.5;
% end

r = zeros(N,1);
rq = 0.5 * ones(N,1); % q trajectory to follow

y_sim(1:p,1) = y_d(1:p,T + 1);
y_sim(2,1) = SOC0;
x(1, 1) = y_d(1,T + 1); % Initial states
x(2, 1) = x_d(2,T + 2); % Initial states
x(3, 1) = SOC0;

% Varying generation / load balance at the node

true_y_sim = y_sim; % true system dynamics without measurement noise, for plots

% Known model: [x_k(1) ... x_k(N+1)] = M [x_ini u(1) ... u(N-1)]^T
% Hard coded; matrix should be
% [I 0 .... 0
% [A B 0 .... 0
% [A^2 AB B 0 .... 0
% etc
for I = 1:N
    for J = 1:N
        if I >= J
            if J == 1
                M(I,J) = 1;
            else
                M(I,J) = -ts*tq;
            end
        end
    end
end

% cvx_solver sedumi

%% Main loop
cost_u(lengthSim) = 0;
cost_y(lengthSim) = 0;
u_sim(2,lengthSim) = 0;
comp_cost(lengthSim) = 0;
for k = 1:lengthSim

    k

    predictedGenLoad  = u_genload(k:k+N-1)';

    % Solve optimisation problem using cvx
    tic
   
    cvx_begin quiet
    cvx_precision best

    variable g(G,1)
    variable u(N*m,1);
    variable y(N*p,1);
    variable sigma_y(Tini*p,1);
    variable q(N,1);

    minimize quad_form(q-rq,Q2expanded) + quad_form(y-r,Qexpanded) + quad_form(u,Rexpanded) + lg * quad_form(g,eye(G)) + ly * norm(sigma_y,1) % Complete cost function

    [Up; Yp; Uf; Yf] * g == [uini; (yini + sigma_y); u; y] % trajectories have to be consistent with initial data

    % constraints on u and y
    for i = 1:m:m*N
        u(i) <= Imax
        u(i) >= -Imax

        u(i+1) == predictedGenLoad(1 + floor(i/2));
    end

    for i = 1:N
        y(i) <= dV
        y(i) >= -dV

        q(i) >= 0
        q(i) <= 1
    end

    if k > 1
        q == M * [y_sim(2,k-1); u(1:m:N*m-m)];
    else
        q == M * [x(3,k); u(1:m:N*m-m)];
    end

    %     q(1) == x(2,k);
    %     for i = 2:N
    %         q(i) == q(i-1) - ts * tq * u(i-1);
    %     end

    cvx_end

    comp_cost(k) = toc;

    if ~isnan(u(1))
        u_sim(1,k) = u(1);
    else
        fprintf('Optimization problem failed to solve \n')
        u_sim(1,k) = 0;
    end
    u_sim(2,k) = u(2);

    % Update system
    x(:, k+1) = A * x(:, k) + B * u_sim(:, k);
    y_sim(:, k) = C * x(:, k) + D * u_sim(:, k) + noiseM*randn(p,1);
    true_y_sim(:, k) = C * x(:, k) + D * u_sim(:, k);

    % Update data vectors - new data at the end of the vector
    uini(1:m,:) = [];
    yini(1:p,:) = [];
    uini = [uini; u_sim(:, k)];
    yini = [yini; y_sim(1, k)];

    % Record costs
    cost_y(k) = quad_form(true_y_sim(1, k)-r(1),Q(1,1)) + quad_form(true_y_sim(2, k)-rq(1),Q(2,2));
    cost_u(k) = R * u_sim(1,k) ^ 2;
end

%% Plots
% figure(1)
% subplot(4,1,1)
% plot(1:k,u_sim(1,1:k))
% title('BESS current')
% 
% subplot(4,1,2)
% plot(1:k,u_sim(2,1:k))
% title('Generation/load balance')
% 
% subplot(4,1,3)
% plot(1:k,true_y_sim(1,1:k))
% title('Node voltage deviation')
% 
% subplot(4,1,4)
% plot(1:k,true_y_sim(2,1:k))
% title('BESS SoC')

figure(1)
subplot(3,1,1)
plot(1:k,u_sim(1,1:k))
title('BESS current')

subplot(3,1,2)
plot(1:k,true_y_sim(1,1:k))
title('Node voltage deviation')

subplot(3,1,3)
plot(1:k,true_y_sim(2,1:k))
title('BESS SoC')

save hdeepc_sq_reg