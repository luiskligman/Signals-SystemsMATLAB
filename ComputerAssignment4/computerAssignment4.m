close all;
clear;
clc;


% PROBLEM 1
n = 30; % Change this to any n

% Approximate EGF: f(t) = sum F_k / k! * t^k
% Since f^(n)(0) = F_n, we simulate it directly

% We'll hardcode Fibonacci values up to n
F = zeros(1, n+1);
F(2) = 1;
for i = 3:n+1
    F(i) = F(i-1) + F(i-2);
end

% f(t) = sum_{k=0}^{n} (F_k / k!) * t^k
% The coefficient of t^n is F_n / n!
% So: F_n = coeff * n!

coeff = F(n+1) / factorial(n);
Fn = coeff * factorial(n);

fprintf('F_%d = %d\n', n, Fn);


% PROBLEM 2
% Fibonacci numbers F_0 to F_15
F = zeros(1, 16);
F(2) = 1;
for i = 3:16
    F(i) = F(i-1) + F(i-2);
end

% Time point to test (can use any small value near 0)
t = 0.1;

% Compute f(t)
f = 0;
for k = 0:15
    f = f + (F(k+1) / factorial(k)) * t^k;
end

% Numerical derivatives using finite difference
h = 1e-5;

% f'(t)
f_plus = 0;
f_minus = 0;
for k = 0:15
    f_plus = f_plus + (F(k+1) / factorial(k)) * (t + h)^k;
    f_minus = f_minus + (F(k+1) / factorial(k)) * (t - h)^k;
end
df = (f_plus - f_minus) / (2*h);

% f''(t)
f_plus2 = 0;
f_minus2 = 0;
for k = 0:15
    f_plus2 = f_plus2 + (F(k+1) / factorial(k)) * (t + h)^k;
    f_minus2 = f_minus2 + (F(k+1) / factorial(k)) * (t - h)^k;
end
d2f = (f_plus2 - 2*f + f_minus2) / (h^2);

% Check identity
check = abs(d2f - (df + f)) < 1e-6;

fprintf("f''(t) ≈ f'(t) + f(t) at t = %.4f: %s\n", t, string(check));


% PROBLEM 3
% Define s range (for possible later evaluation or plotting)
s = linspace(-5, 5, 1000);  % Just a dummy range if needed

% Manually compute F(s) from:
% (s^2 - s - 1) * F(s) = 1  -->  F(s) = 1 / (s^2 - s - 1)

% Example: Evaluate F(s) at a few points for verification
s_vals = [2, 3, 4];
for i = 1:length(s_vals)
    s_i = s_vals(i);
    Fs = 1 / (s_i^2 - s_i - 1);
    fprintf('F(%d) = %.5f\n', s_i, Fs);
end


% PROBLEM 4
% Define numerator and denominator
num = [1];           % Numerator: constant → no finite zeros
den = [1 -1 -1];     % Denominator: s^2 - s - 1

% Compute zeros and poles
zeros_F = roots(num);   % Zeros of F(s)
poles_F = roots(den);   % Poles of F(s)

% Display results
fprintf('Zeros of F(s):\n');
if isempty(zeros_F)
    fprintf('  None (numerator is constant)\n');
else
    disp(zeros_F);
end

fprintf('Poles of F(s):\n');
disp(poles_F);


% PROBLEM 5a
% Approximate f(t) using EGF with F_0 to F_5
% Define Fibonacci numbers F_0 to F_5
F = [0 1 1 2 3 5];

% Time range
t = linspace(-2, 2, 400);
f_approx = zeros(size(t));

% Compute EGF: f(t) = sum_{k=0}^5 (F_k / k!) * t^k
for k = 0:5
    f_approx = f_approx + (F(k+1) / factorial(k)) * t.^k;
end

% Plot the result
figure;
plot(t, f_approx, 'b', 'LineWidth', 2);
xlabel('t');
ylabel('f(t)');
title('Problem 5a: EGF Approximation using F_0 to F_5');
grid on;

% PROBLEM 5b
% Plot exact f(t) from inverse Laplace
phi1 = (1 + sqrt(5)) / 2;
phi2 = (1 - sqrt(5)) / 2;

% Unit step: u(t) = 1 when t >= 0
u = double(t >= 0);

% Compute f(t)
f_exact = (exp(phi1 * t) - exp(phi2 * t)) / sqrt(5) .* u;

% Plot exact f(t)
figure;
plot(t, f_exact, 'r', 'LineWidth', 2);
xlabel('t');
ylabel('f(t)');
title('Problem 5b: Exact f(t) from Inverse Laplace Transform');
grid on;


% PROBLEM 6
n = 12;  % or any n ≥ 0

phi1 = (1 + sqrt(5)) / 2;
phi2 = (1 - sqrt(5)) / 2;

F_n = (phi1^n - phi2^n) / sqrt(5);

fprintf('F_%d = %.0f\n', n, F_n);


% PROBLEM 7
n = 19;

phi1 = (1 + sqrt(5)) / 2;
phi2 = (1 - sqrt(5)) / 2;

F_19 = (phi1^n - phi2^n) / sqrt(5);

fprintf('F_%d = %.0f\n', n, F_19);