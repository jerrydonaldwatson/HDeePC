%% Example 1c (sysID + MPC)
% Jeremy Watson, 2025
%
% This code accompanies our paper: "Hybrid Data-Enabled Predictive
% Control: Incorporating model knowledge into the DeePC"
%
% Contact: jeremy.watson@canterbury.ac.nz in case of any queries.

clear
close all
clc

rng(1)

%% System Set Up
load config

%% System ID
m = 2;          % Number of inputs
p = 2;          % Number of outputs
n = 3;          % Number of states

% Initial simulation to gather data
u = (rand(m, T+1) - 0.5);
y = zeros(p, T);
x(:,1) = rand(n,1);
for k = 1:T+1
    if u(1,k) >= 0
        x(:, k+1) = [a * sin(x(1:2,k)) + b * (x(1:2,k).*[u(:,k) + 1]); x(3,k)] + B1*u(:, k);
    else
        x(:, k+1) = [a * sin(x(1:2,k)) + b * (x(1:2,k).*[u(:,k) + 1]); x(3,k)] + B2*u(:, k);
    end
    y(1:p, k) = C(1:p, :)*x(:,k) + noiseM*randn(p,1);
end

sysID = n4sid(iddata(y', u', ts), 3);

R = diag([R R]);

y_sim(1:p,1) = y(1:p,T + 1);
x_d = x;
clear x
x(:, 1) = x_d(:,T + 2); % Initial states

true_y_sim = y_sim; % true system dynamics without measurement noise, for plots


%% Main loop
cost_u(lengthSim) = 0;
cost_y(lengthSim) = 0;
tic
for k = 1:lengthSim

    k

    % Solve optimisation problem using cvx
    cvx_begin quiet
    cvx_precision best

    variable u(m,N)
    variable x_k(n,N+1)
    variable y(p,N)

    if useSlack
        variable sigma_x(n,N)
        variable sigma_y(p,N)
    end

    cost = 0;
    for i = 1:N
        cost = cost + (y(:,i))'*Q*(y(:,i)) + u(:,i)'*R*u(:,i);
    end

    if useSlack
        for i = 1:N
            cost = cost + lx*sigma_x(:,i)'*sigma_x(:,i) + ly*sigma_y(:,i)'*sigma_y(:,i);
        end
    end

    minimize cost

    y(:,1) == C*x(:,k);

    for i = 1:N
        if useSlack
            x_k(:,i+1) == sysID.A * x_k(:,i) + sysID.B * u(:,i) + sigma_x(:,i);
            y(:,i) == sysID.C * x_k(:,i) + sysID.D * u(:,i) + sigma_y(:,i);
        else
            x_k(:,i+1) == sysID.A * x_k(:,i) + sysID.B * u(:,i);
            y(:,i) == sysID.C * x_k(:,i) + sysID.D * u(:,i);
        end
    end

    for i = 1:N
        u(:,i) <= Umax
        u(:,i) >= -Umax
    end

    for i = 1:N
        y(1,i) <= dX1
        y(1,i) >= -dX1
        y(2,i) <= dX3
        y(2,i) >= -dX3
    end

    cvx_end

    % Assign outputs for next time-step (at present, recalculating every step...)
    if ~isnan(u(1,1)) && ~(cvx_status == "Failed" || cvx_status == "Infeasible")
        u_sim(:,k) = u(1:m,1);
    else
        fprintf('Optimization problem failed to solve \n')
        pause
        u_sim(:,k) = u(1:m,1);
    end


    % Update system
    if u_sim(1,k) >= 0
        x(:, k+1) = [a * sin(x(1:2,k)) + b * (x(1:2,k).*[u_sim(:,k) + 1]); x(3,k)] + B1 * u_sim(:, k);
    else
        x(:, k+1) = [a * sin(x(1:2,k)) + b * (x(1:2,k).*[u_sim(:,k) + 1]); x(3,k)] + B2 * u_sim(:, k);
    end
    true_y_sim(:, k) = C * x(:, k) + D * u_sim(:, k);

    % Record costs
    cost_y(k) = quad_form(true_y_sim(:, k),Q);
    cost_u(k) = R(1,1) * u_sim(1,k) ^ 2 + R(2,2) * u_sim(2,k) ^ 2;
end
toc


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


save bess_nl_sid_mpc