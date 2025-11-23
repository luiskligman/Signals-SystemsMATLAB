% Computer Assignment 8 --- Luis Kligman
clear all;
close all;
clc;

%% Problem 1a: Plot v_s(t) using 100 harmonics
T = 2e-3;  % Period in seconds
w0 = 2*pi/T;  % 1000*pi rad/s
K = 100;  % Number of haromincs on each side
k = -K:K;

% Coefficient vector
ck = zeros(size(k));

% Fill ck
for i = 1:length(k)
    ki = k(i);
    if ki == 0
        ck(i) = 5/2;  % c0
    else
        ck(i) = (5./(pi*ki)) * sin(ki*pi/2);
    end
end

% Time axis for plotting several periods
t = linspace(-5e-3, 5e-3, 4000);

% Reconstruct v_s(t) using complex exponential series
vs = zeros(size(t));
for i = 1:length(k)
    vs = vs + ck(i) .* exp(1j * k(i) * w0 .* t);
end

figure;
plot(t*1e3, real(vs), 'LineWidth', 1.2);
grid on;
xlabel('t (ms)');
ylabel('v_s(t) (V)');
title('Reconstruction of v_s(t) with 100 harmonics');

%% Problem 1b: Magnitude Spectrum of v_s(t)

% Frequency range 
w_min = -10000*pi;
w_max = 10000*pi;

% Harmonics that fall in this range
k_spec = -10:10;
w_spec = k_spec * w0;

% Compute |c_k| for these k values:
ck_spec = zeros(size(k_spec));
for i = 1:length(k_spec)
    ki = k_spec(i);
    if ki == 0
        ck_spec(i) = 5/2;
    else
        ck_spec(i) = (5/(pi*ki)) * sin(ki*pi/2);
    end
end

% Plot the magnitude spectrum
figure;
stem(w_spec, abs(ck_spec), 'filled');
grid on;
xlabel('\omega (rad/s)');
ylabel('|c_k|');
title('Magnitude Spectrum of v_s(t)');
xlim([w_min, w_max]); 

%% Problem 2a: Plot v_out(t) using 100 harmonics

% Solved transfer function on paper; H(jw)
H = @(w) 4000 ./ (8000 + 1j*w);

% Output Fourier coefficients d_k
dk = ck .* H(k * w0);

% Reconstruct v_out(t)
vout = zeros(size(t));
for i = 1:length(k)
    vout = vout + dk(i) .* exp(1j * k(i) *w0 .* t);
end

figure;
plot(t*1e3, real(vout), 'LineWidth', 1.8);
grid on;
xlabel('t (ms)');
ylabel('v_{out}(t) (V)');
title('Reconstruction of v_{out}(t) with 100 harmonics');

%% Problem 2b: Magnitude spectrum of v_out(t) 
k_spec = -10:10;
w_spec = k_spec * w0;

dk_spec = zeros(size(k_spec));
for i = 1:length(k_spec)
    ki = k_spec(i);
    if ki == 0
        dk_spec(i) = H(0) * (5/2);
    else
        ck_i = (5/(pi*ki)) * sin(ki*pi/2);
        dk_spec(i) = H(ki*w0) * ck_i;
    end
end

figure;
stem(w_spec, abs(dk_spec), 'filled');
grid on;
xlabel('\omega (rad/s)');
ylabel('|d_k|');
title('Magnitude Spectrum of v_{out}(t)');
xlim([-10000*pi 10000*pi]);
