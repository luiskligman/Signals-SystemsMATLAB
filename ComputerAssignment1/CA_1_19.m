% Time vector
t = -5:0.01:5;  % from -5s to 5s with step size 0.01

% Define the unit step function
u = @(x) double(x >= 0);

% Define signals
x1 = -6 * u(t+3);
x2 = 10 * u(t - 4);
x3 = 4 * u(t + 2) - 4 * u(t - 2);

% Plot x1(t)
figure;
plot(t, x1, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x1(t)');
title('x1(t) = -6u(t+3)');
ylim([-7 2]);

% Plot x2(t)
figure;
plot(t, x2, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x2(t)')
title('x2(t) = 10u(t-4)');
ylim([-1 11]);

% Plot x3(t)
figure;
plot(t, x3, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x3(t)');
title('x3(t) = 4u(t+2) - 4u(t-2)');
ylim([-1 5]);
