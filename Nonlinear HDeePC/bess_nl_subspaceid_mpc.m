%% BESS Simulation (DeePC)
clear
close all
clc


rng(1)

%% System Set Up
load config

% Define System Matrices

A = [1-ts*(1/Rload/Cdc), ts/Cdc 0; % voltage
    -ts*(1/Lline) 1-ts*Rline/Lline 0;
    0, 0, 1]; % SoC

B1 = [ts/Cdc ts/Cdc;
    0 0;
    -ts*tq/e 0];

B2 = [ts/Cdc ts/Cdc;
    0 0;
    -ts*tq*e 0];

C = [1, 0, 0;
    0 0 1];

D = [0 0;
    0 0];

%% System ID
m = 2;          % Number of inputs
p = 2;          % Number of outputs
n = 3;          % Number of states

% Initial simulation to gather data
u = 10 .* (rand(m, T+1) - 0.5);
y = zeros(p, T);
x(1,1) = rand(1,1);
x(2,1) = 0;
x(3,1) = rand(1);
SOC0 = x(3,1)
for k = 1:T+1
    
    if u(1,k) >= 0
        x(:, k+1) = A*x(:, k) + B1*u(:, k);
    else
        x(:, k+1) = A*x(:, k) + B2*u(:, k);
    end

    y(:, k) = C*x(:,k) + noiseM*randn(p,1);
end

sysID = n4sid(iddata(y', u', ts), 3);

% Trajectory to try to follow
r(1,1) = 0;
r(2,1) = 0.5;

R = diag([R 0]);

y_sim(1:p,1) = y(1:p,T + 1);
y_sim(2,1) = SOC0;
x(1, 1) = x(1,T + 2); % Initial states
x(2, 1) = x(2,T + 2); % Initial states
x(3, 1) = x(3,T + 2);

x
% Varying generation / load balance at the node

true_y_sim = y_sim; % true system dynamics without measurement noise, for plots

%% Main loop
cost_u(lengthSim) = 0;
cost_y(lengthSim) = 0;
tic
for k = 1:lengthSim

    k

    if k >= lengthSim/2
        r = 0.72;
    else
        r = 0.7;
    end

    predictedGenLoad  = u_genload(k:k+N-1)';

    % Solve optimisation problem using cvx
    cvx_begin quiet
    cvx_precision best

    variable u(m,N)
    variable x_k(n,N+1)
    variable y(p,N)

    cost = 0;
    for i = 1:N
        cost = cost + (y(:,i)-r)'*Q*(y(:,i)-r) + u(:,i)'*R*u(:,i);
    end

    minimize cost

    y(:,1) == C*x(:,k);

    for i = 1:N
        x_k(:,i+1) == sysID.A * x_k(:,i) + sysID.B * u(:,i);
        y(:,i) == sysID.C * x_k(:,i) + sysID.D * u(:,i);
    end

    u(2,:) == predictedGenLoad';

    for i = 1:N
        u(1,i) <= Imax
        u(1,i) >= -Imax
    end

    for i = 1:N
        y(1,i) <= dV
        y(1,i) >= -dV
    end

    cvx_end

    % Assign outputs for next time-step (at present, recalculating every step...)
    if ~isnan(u(1,1))
        u_sim(:,k) = u(1:m,1);
    else
        fprintf('Optimization problem failed to solve \n')
        u_sim(:,k) = 0;
    end


    % Update system
    if u_sim(1,k) >= 0
        x(:, k+1) = A * x(:, k) + B1 * u_sim(:, k);
    else
        x(:, k+1) = A * x(:, k) + B2 * u_sim(:, k);
    end
    y_sim(:, k) = C * x(:, k) + D * u_sim(:, k) + noiseM*2*(rand(p,1) - 0.5);
    true_y_sim(:, k) = C * x(:, k) + D * u_sim(:, k);

    % Record costs
    cost_y(k) = quad_form(true_y_sim(:, k)-r,Q);
    cost_u(k) = R(1,1) * u_sim(1,k) ^ 2;
end
toc


% Plot input, predicted output, and reference
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
hold on
plot(1:k, [0.70*ones(floor(k/2),1); 0.72*ones(floor(k/2),1)]) % plot reference

save bess_nl_sid_mpc