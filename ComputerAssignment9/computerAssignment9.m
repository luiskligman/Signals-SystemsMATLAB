% Luis Kligman --- Computer Assignment 9
clear all;
close all;
clc;

%% Parameters
A = 5e-3;  % 5 mA
T = 3;  % Period (seconds)
R1 = 500;
R2 = 2000;
C = (1/3) * 1e-3;  % 1/3 mF

tau = C*(R1 + R2);  % time constant = 5/6 seconds

% Define I_s(omega)
I_s = @(w) ...
    (w==0).*(A*T/2) + ...
    (w~=0).* ( (8*A ./ (T*w.^2)) .* sin(T*w/4).^2 );

% Transfer function H(omega)
H = @(w) R1 ./ (1 + 1j*w*tau);

% Output spectrum V_out(omega)
V_out = @(w) H(w).*I_s(w);

%% Part D: plot functions with a numeric integration from -inf to inf

% Time vector
t = linspace(-2,5,400);

i_s = zeros(size(t));
v_out = zeros(size(t));

wMax = 400;

for k = 1:length(t)
    ti = t(k);
    i_s(k) = real(integral(@(w) I_s(w) .* exp(1j*w*ti), -wMax, wMax) / (2 * pi));
    v_out(k) = real(integral(@(w) V_out(w) .* exp(1j*w*ti), -wMax, wMax) / (2 * pi));
end

% Plot i_s(t)
figure;
plot(t, i_s, 'LineWidth', 2);
xlabel('t');
ylabel('i_s(t)');
title('Reconstructed Input Current');

% Plot v_out(t)
figure;
plot(t, v_out, 'LineWidth', 2);
xlabel('t');
ylabel('v_{out}(t)');
title('Reconstructued Output Voltage');

%% Part E: Plot one-sided energy spectral density of v_out(t)

w = linspace(0, wMax, 2000);  % rad/s

% Two-sided energy spectral density
Sv_two = abs(V_out(w)).^2;

% One-sided energy spectral density
% double everything except DC 
Sv_one = 2 * Sv_two;
Sv_one(1) = Sv_two(1);  % Do not double the DC component

% Plot one-sided energy spectral density
figure;
plot(w, Sv_one, 'LineWidth', 2);
xlabel('\omega (rad/s)');
ylabel('S_v^{(1)}(\omega)');
title('One-Sided Energy Spectral Density of v_{out}(t)');
grid on;

