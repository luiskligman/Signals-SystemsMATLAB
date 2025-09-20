% ELCT 222 - Computer Assignment 3

clear;
clc;
close all;

% Constants
RC_values = [0.1, 1];  % Case (i) and (ii)
f = linspace(0, 10, 1000);  % Frequency range (Hz)
t = linspace(0, 2*pi, 1000);  % Time vector for x(t) and y(t)

% Input signal
x = 1 + cos(5*t) + cos(10*t) + cos(15*t);

% Part 1a and 1b: Frequency response plots
for idx = 1:2
    RC = RC_values(idx);

    % Frequency response H(f)
    H_f = 1 ./ (1 + 1j*2*pi*f*RC);
    mag_H = abs(H_f);
    phase_H = angle(H_f);

    % Plot magnitude
    figure;
    subplot(2, 1, 1);
    plot(f, mag_H, 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('|H(f)|');
    title(['Magnitude Response for RC = ', num2str(RC)]);
    grid on;

    % Plot phase
    subplot(2,1,2);
    plot(f, phase_H, 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Phase(H(f)) [rad]');
    title(['Phase Response for RC = ', num2str(RC)]);
    grid on;
end

% Part 1c: Compute output y(t) using frequency response at 0, 5, 10, 15 Hz
frequencies = [0, 5, 10, 15];  % Frequencies at which to compute output

for idx = 1:2
    RC = RC_values(idx);

    % Compute frequency response values
    H = 1 ./ (1 + 1j*2*pi*frequencies*RC);

    % Output signal y(t)
    y = real(H(1)) + ...
        abs(H(2))*cos(5*t + angle(H(2))) + ...
        abs(H(3))*cos(10*t + angle(H(3))) + ...
        abs(H(4))*cos(15*t + angle(H(4)));

    % Store for later plot
    if idx == 1
        y_RC_01 = y;
    else
        y_RC_1 = y;
    end
end

% Part 1d: Plot |x(t)| and |y(t)| for both RC values
for idx = 1:2
    RC = RC_values(idx);

    if idx == 1
        y = y_RC_01;
    else
        y = y_RC_1;
    end

    figure;
    plot(t, abs(x), 'b', 'LineWidth', 1.5);
    hold on;
    plot(t, abs(y), 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('|x(t)|', '|y(t)|');
    title(['Input vs Output Magnitude for RC = ', num2str(RC)]);
    grid on;
end
