%% Example 1c, DeePC
% Jeremy Watson, 2025
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.

clear all
close all
clc

rng(1)

%% System Set Up
load config

%% DeePC Controller
m = 2;          % Number of inputs
p = 2;          % Number of outputs
n = 3;          % Number of states
for I = 1:N
    Qexpanded(2*I-1,2*I-1) = Q(1,1);
    Qexpanded(2*I,2*I) = Q(2,2);
end

Rexpanded = diag(R * ones(1,N*m));

% Initial simulation to gather data
u_d = (rand(m, T+1) - 0.5);
y_d = zeros(p, T);
x_d(:,1) = rand(n,1);
for k = 1:T+1
    if u_d(1,k) >= 0
        x_d(:, k+1) = [a * sin(x_d(1:2,k)) + b * (x_d(1:2,k).*[u_d(:,k) + 1]); x_d(3,k)] + B1*u_d(:, k);
    else
        x_d(:, k+1) = [a * sin(x_d(1:2,k)) + b * (x_d(1:2,k).*[u_d(:,k) + 1]); x_d(3,k)] + B2*u_d(:, k);
    end
    y_d(1:p, k) = C(1:p, :)*x_d(:,k) + noiseM*randn(p,1);
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

y_sim(1:p,1) = y_d(1:p,T + 1);
x(:, 1) = x_d(:,T + 2); % Initial states

true_y_sim = y_sim; % true system dynamics without measurement noise, for plots

% cvx_solver sedumi

%% Main loop
cost_u(lengthSim) = 0;
cost_y(lengthSim) = 0;
u_sim(2,lengthSim) = 0;
tic
for k = 1:lengthSim

    k

    % Solve optimisation problem using cvx
    cvx_begin quiet
    cvx_precision best

    variable g(G,1)
    variable u(N*m,1);
    variable y(N*p,1);
    variable sigma_y(Tini*p,1);

    minimize (y'*Qexpanded*y + u'*Rexpanded*u + lg * norm(g,1) + ly * norm(sigma_y,1)) % Complete cost function

    [Up; Yp; Uf; Yf] * g == [uini; (yini + sigma_y); u; y] % trajectories have to be consistent with initial data

    % constraints on u and y
    for i = 1:m:m*N
        u(i) <= Umax
        u(i) >= -Umax

        u(i+1) <= Umax
        u(i+1) >= -Umax
    end
    for i = 1:p:p*N
        y(i) <= dX1
        y(i) >= -dX1
        y(i+1) <= dX3
        y(i+1) >= -dX3
    end

    cvx_end

    % Assign outputs for next time-step (at present, recalculating every step...)
    u_sim(1,k) = u(1);
    u_sim(2,k) = u(2);
    if ~isnan(u(1))  && ~(cvx_status == "Failed" || cvx_status == "Infeasible")
    else
        fprintf('Optimization problem failed to solve \n')

        pause
    end


    % Update system (nonlinear system depending on the sign of u_1)
    if u_sim(1,k) >= 0
        x(:, k+1) = [a * sin(x(1:2,k)) + b * (x(1:2,k).*[u_sim(:,k) + 1]); x(3,k)] + B1 * u_sim(:, k);
    else
        x(:, k+1) = [a * sin(x(1:2,k)) + b * (x(1:2,k).*[u_sim(:,k) + 1]); x(3,k)] + B2 * u_sim(:, k);
    end
    y_sim(:, k) = C * x(:, k) + D * u_sim(:, k) + noiseM*randn(p,1);
    true_y_sim(:, k) = C * x(:, k) + D * u_sim(:, k);

    % Update data vectors - new data at the end of the vector
    uini(1:m,:) = [];
    yini(1:p,:) = [];
    uini = [uini; u_sim(:, k)];
    yini = [yini; y_sim(:, k)];

    % Record costs
    cost_y(k) = quad_form(true_y_sim(:, k),Q);
    cost_u(k) = R * u_sim(1,k) ^ 2 + R * u_sim(2,k) ^ 2;
end
toc

%% Plots
figure(1)
subplot(3,1,1)
plot(1:k,u_sim(:,1:k))
legend('u_1, u_2')
title('Inputs')

subplot(3,1,2)
plot(1:k,true_y_sim(1,1:k))
legend('y_1')

subplot(3,1,3)
plot(1:k,true_y_sim(2,1:k))
legend('y_2')


save deepc