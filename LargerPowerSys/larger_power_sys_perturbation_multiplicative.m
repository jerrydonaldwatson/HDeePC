%% Larger (power) system example (Example 3)
% Jeremy Watson, 2025
%
% Acknowledgement:
% power_grid_uXX.mat and the example setting come from
% U. Wasekar, J. Watson "Monte-Carlo analysis of interlinking converter
% modelling and control in hybrid AC/DC networks", PMAPS2024.
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
% Example 3 (Perturbation resulting in model mismatch)
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.

clear all
close all
clc

%% System Set Up
lengthSim = 10; % Timestep is 0.1 s
noiseM = 1e-5; % Uniform distribution with this limit
perturbationLevel = 0.2; % = standard deviation of perturbationLevel * 100 %

% Load System Matrices
load power_grid_u00.mat
Atrue = power_grid_ssd.A;
Btrue = power_grid_ssd.B;
Btrue(:,3:6) = Btrue(:,3:6) .* 1e3; % Scale units to kW
Ctrue = eye(length(Atrue));
Dtrue = zeros(length(Atrue), size(Btrue,2));

% Perturb system matrices multiplicatively (structure preserving)
rng(1)
A = power_grid_ssd.A .* (1 + perturbationLevel*randn(size(Atrue)));
B = power_grid_ssd.B .* (1 + perturbationLevel*randn(size(Btrue)));
B(:,3:6) = B(:,3:6) .* 1e3; % Scale units to kW
C = eye(length(A));
D = zeros(length(A), size(B,2));

% Loop adjusting number of known states
for p_u = [0 9 17 34 51 60 68]

    p_u

    clearvars -except A B C D Atrue Btrue Ctrue Dtrue lengthSim noiseM p_u

    rng(0)

    %% HDeePC Controller

    % Adjust the number of unknown states n_u:
    % n_u = 0 => MPC
    % 0 < n_u < n => HDeePC
    % n_u = n => DeePC
    % Also ensure that the p_k is set correctly

    m = size(B,2);          % Number of inputs
    p = size(C,1);          % Number of outputs
    % p_u = 20;               % Number of outputs for which the equation is known
    p_k = p - p_u;          % Number of outputs for which the equation is unknown
    n = size(A,1);          % Number of states
    n_u = p_u;              % Number of unknown states (since C = I)
    n_k = n - n_u;          % Number of known states

    % Parameters
    Tini = 68;      % Past data for learning the system
    N = 20;         % Horizon
    T = 500;        % Data length;
    reducedTol = 1e-3; % Reduced tolerance for large-scale problem

    % Cost function weights
    Q = 0.001; % All states are currently weighted the same
    R = 0.05; 
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
    A_u = A(1:n_u,1:n_u);
    A_f = A(1:n_u,n_u+1:end);
    A_c = A(n_u+1:end,1:n_u);
    A_k = A(n_u+1:end,n_u+1:end);
    B_c = B(1:n_u,:);
    B_k = B(n_u+1:end,:);
    C_u = C(1:p_u,1:n_u);
    C_f = C(1:p_u,n_u+1:end);
    C_c = C(p_u+1:end,1:n_u);
    C_k = C(p_u+1:end,n_u+1:end);
    D_k = D(p_u+1:end,:);

    A_y = A_c(1:n_k,:); % remember to adjust if y_u =/= x_u
    C_y = C_c(1:p_k,:);

    %% Initial simulation to gather data
    u_d = (rand(m, T+1) - 0.5);
    y_d = zeros(p_u, T); % p_u outputs for data-based part of model
    y_d_all = zeros(p, T); % used to record all outputs for info
    x_d(:,1) = rand(n,1);
    for k = 1:T+1
        x_d(:, k+1) = Atrue*x_d(:, k) + Btrue*u_d(:, k);

        noise = noiseM*2*(rand(p,1) - 0.5);
        y_d(:, k) = Ctrue(1:p_u,:)*x_d(:,k) + noise(1:p_u,1);
        y_d_all(:, k) = Ctrue*x_d(:,k) + noise;
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
    x(:, 1) = (rand(n,1) - 0.5) ./ 1e2; % start system off perturbed

    cvx_solver mosek

    %% Main loop
    cost_u(lengthSim) = 0;
    cost_y(lengthSim) = 0;
    tic
    for k = 1:lengthSim

        % Solve optimisation problem using cvx
        cvx_begin quiet % quiet for timing

        variable g(G,1)
        variable u(N*m,1);
        variable y_u(N*p_u,1);
        variable y_k(N*p_k,1);
        variable x_k(n_k,N+1);
        variable sigma_y(Tini*p_u,1);
        variable sigma_u(Tini*m,1);

        expr_yk = reshape(y_k - r_k, [], 1); % This reformulation is to prevent a vec() overload error on some installations
        expr_yu = reshape(y_u - r, [], 1);

        minimize( sum_square(chol(Qexpanded) * expr_yu) ...
        + sum_square(chol(Q * eye(N*p_k)) * expr_yk) ...
        + sum_square(chol(Rexpanded) * u) ...
        + lg * norm(g,1) ...
        + lu * norm(sigma_u,2) ...
        + ly * norm(sigma_y,2) )

        % Data-based part
        if n_u > 0
            [Up; Yp; Uf; Yf] * g == [(uini + sigma_u); (yini + sigma_y); u; y_u] % trajectories have to be consistent with initial data
        end

        % Constraints on u and (possibly) y
        for i = 1:N*m
            u(i) <= 1
            u(i) >= -1
        end

        % MPC part
        x_k(:,1) == x(n_u+1:end, k);
        for i = 1:N
            x_k(:,i+1) == A_y * y_u((1:p_u)+(i-1)*p_u) + A_k * x_k(:,i) + B_k * u((1:m)+(i-1)*m); % note that y_u(1:2) = x_u(1:2)
            y_k((1:p_k)+(i-1)*p_k) == C_y * y_u((1:p_u)+(i-1)*p_u) + C_k * x_k(:,i) + D_k * u((1:m)+(i-1)*m);
        end

        cvx_end

        % Assign outputs for next time-step (at present, recalculating every step...)
        if ~isnan(u(1)) && cvx_status ~= "Infeasible"  && (cvx_status ~= "Failed" || cvx_slvtol < reducedTol)
            u_sim(:,k) = u(1:m,1);
        else
            fprintf('Optimization problem failed to solve \n')
            cvx_slvtol
            u_sim(:,k) = u(1:m,1);
        end

        % Update system
        x(:, k+1) = Atrue * x(:, k) + Btrue * u_sim(:, k);
        y_sim(:, k) = Ctrue * x(:, k) + Dtrue * u_sim(:, k) + noiseM*2*(rand(p,1) - 0.5);
        true_y_sim(:, k) = Ctrue * x(:, k) + Dtrue * u_sim(:, k);

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

end