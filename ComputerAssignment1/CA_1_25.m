close all; 
clear;
clc;

u = @(x) double(x >= 0);  % Unit step

% (a) --- x1(t) = 100e^{-2t} u(t)  -> tau = 1/2
t = 0:0.001:3;
x1 = 100 * exp(-2 * t).*u (t);
figure;
plot(t,x1,'LineWidth',2);
grid on;
xlabel('t (s)');
ylabel('x_1(t)');
title('x_1(t) = 100e^{-2t} u(t)');

% (b) --- x2(t) = -10e^{-0.1t} u(t)  --> tau = 10
t = 0:0.1:70;  % ~7 * tau to show decay to near zero
x2 = -10 * exp(-0.1 * t).*u(t);
figure;
plot(t, x2, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_2(t)');
title('x_2(t) = -10e^{-0.1} u(t)');

% (c) --- x3(t) = -10e^{-0.1t} u(t-5) --> starts at t = 5
t = 0:0.1:70;
x3 = -10 * exp(-0.1 * t).*u(t-5);
figure;
plot(t, x3, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_3(t)');
title('x_3(t) = -10e^{-0.1t} u(t-5)')

% (d) --- x4(t) = 10(1 - e^{-10^3 t}) u(t) --> tau = 1/1000
t = 0:1e-5:0.01;  % very fast rise; show first 10 ms
x4 = 10 * (1 - exp(-1e3*t)).*u(t);
figure;
plot(t, x4, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_4(t)');
title('x_4(t) = 10(1-e^{-10^3 t} u(t)')

% (e) --- x5(t) = 10e^{-0.2 (t - 4)} u(t)  --> tau = 5, active from t≥0
t = 0:0.05:40;
x5 = 10*exp(-0.2*(t - 4)).*u(t);
figure;
plot(t, x5, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_5(t)');
title('x_5(t) = 10e^{-0.2(t - 4)} u(t)');

% (f) --- x6(t) = 10e^{-0.2 (t - 4)} u(t - 4)  --> tau = 5, starts at t = 4
t = 0:0.05:40;
x6 = 10*exp(-0.2*(t-4)).*u(t-4);
figure;
plot(t, x6, 'LineWidth', 2);
grid on;
xlabel('t (s)');
ylabel('x_6(t)');
title('x_6(t) = 10e^{-0.2(t-4)} u(t-4)');