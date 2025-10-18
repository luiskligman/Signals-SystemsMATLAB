% Computer Assignment 5 -- Luis Kligman

clear;
clc;
close all;

tspan = [0 10];  % simulation time
titles = {
    'Case 1: R=1/2, L=1, C=1'
    'Case 2: R=1, L=1, C=0'
    'Case 3: R=1, L=0, C=1'
    'Case 4: R=INFINITY, L=1, C=1'
};

cases = {
    struct('R', 0.5, 'L', 1, 'C', 1),   % Case 1
    struct('R', 1,   'L', 1, 'C', 0),   % Case 2 (RL)
    struct('R', 1,   'L', 0, 'C', 1),   % Case 3 (RC)
    struct('R', inf, 'L', 1, 'C', 1)    % Case 4 (LC open circuit)
};

figure; hold on;

for i = 1:4
    R = cases{i}.R;
    L = cases{i}.L;
    C = cases{i}.C;

    if isinf(R)
        % LC: L * d²y/dt² + (1/C)*y = input
        f = @(t, y) [y(2); (1/L)*(1 - y(1)/C)];
        y0 = [0; 0];
    elseif C == 0
        % RL: L * dy'/dt + R*y = 1 (step input)
        f = @(t, y) (1/L)*(1 - R*y);
        [t, y] = ode45(@(t, y) f(t, y), tspan, 0);
        plot(t, y, 'LineWidth', 1.5);
        continue;
    elseif L == 0
        % RC: dy/dt = (1/(R*C))*(1 - y)
        f = @(t, y) (1/(R*C))*(1 - y);
        [t, y] = ode45(@(t, y) f(t, y), tspan, 0);
        plot(t, y, 'LineWidth', 1.5);
        continue;
    else
        % RLC: L*y'' + R*y' + (1/C)*y = 1 (step input)
        f = @(t, y) [y(2); (1/L)*(1 - R*y(2) - y(1)/C)];
        y0 = [0; 0];
    end

    [t, y] = ode45(f, tspan, y0);
    plot(t, y(:,1), 'LineWidth', 1.5);
end

legend(titles, 'Location', 'southeast');
xlabel('Time (s)');
ylabel('Step Response y_{step}(t)');
title('Step Responses for Different RLC Configurations');
grid on;