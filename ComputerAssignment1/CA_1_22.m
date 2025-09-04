% Time vector
t = -4:0.01:4;

% Step, ramp, and rect functions
u = @(x) double(x >= 0);
r = @(x) x .* u(x);
rect = @(x) double(abs(x) <= 0.5);  % rect is 1 if |x| <= 0.5

% Signals
x1 = 5 * r(t + 2) - 5 * r(t);
x2 = 5 * r(t + 2) - 5 * r(t) - 10 * u(t);
x3 = 10 - 5 * r(t + 2) + 5 * r(t);
x4 = 10 * rect((t + 1) / 2) - 10 * rect((t - 3) / 2);
x5 = 5 * rect((t - 1) / 2) - 5 * rect((t - 3) / 2);

% Plotting each function
figure;
plot(t, x1, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_1(t)');
title('x_1(t) = 5r(t+2)-5r(t)');
ylim([-1 11]);

figure;
plot(t, x2, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_2(t)');
title('x_2(t) = 5r(t+2)-5r(t)-10u(t)');
ylim([-1 11]);

figure;
plot(t, x3, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_3(t)');
title('x_3(t) = 10-5r(t+2)+5r(t)');
ylim([-1 11]);

figure;
plot(t, x4, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_4(t)');
title('x_4(t) = 10rect((t+1)/2)-10rect((t-3)/2)');
ylim([-11 11]);

figure;
plot(t, x5, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_5(t)');
title('x_5(t) = 5rect((t-1)/2)-5rect((t-3)/2)');
ylim([-6 6]);