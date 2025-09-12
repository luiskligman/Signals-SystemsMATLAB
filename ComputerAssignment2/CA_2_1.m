% ELCT 222 --- Computer Assignment 2
% Problem 1: Analytical solution of RC circuit response to a rectangular pulse input

% Define a time vector from 0 to 5 seconds with 100 evenly spaced points.
t = linspace(0, 5, 1000);

% Define the input x(t): a rectangular pulse of height 2 from 
% t = 0 to t = 1
% Use logical indexing: (t >= 0 & t < 1) evaluates to 1 in that interval
% and 0 elsewhere
x = 2 * (t >= 0 & t < 1);

% Define the two RC values to analyze:
% Case 1: RC = 0.1 --> faster system
% Case 2: RC = 1   --> slower system
RC_values = [0.1, 1];

% Preallocate a mtrix to store output responses for each RC case
% Each row of y_outputs will hold y(t) for one RC value
y_outputs = zeros(length(RC_values), length(t));

% Loop over each RC value to compute and plot the corresponding output y(t)
for k = 1:length(RC_values)
    RC = RC_values(k);  % Current RC value

    % Initialize output vector y(t) to all zeros for this case
    y = zeros(size(t));

    % --- FIRST INTERVAL: 0 ≤ t < 1 ---
    % During this time, the input x(t) = 2 (constant)
    % Solve: RC * dy/dt + y = 2
    % Homogeneous solution + particular solution gives:
    % y(t) = 2(1 - exp(-t / RC)) with initla condition y(0) = 0
    idx1 = t < 1;  % Indices where t < 1
    y(idx1) = 2 * (1 - exp(-t(idx1) / RC));

    % --- SECOND INTERVAL: t ≥ 1 ---
    % At t = 1, the input x(t) drops to 0.
    % Now solve: RC * dy/dt + y = 0
    % Use y(1) from the previous step as initial codition for continuity
    y1 = 2 * (1 - exp(-1 / RC));  % This is y(t=1)

    % Output continues as exponential decay:
    % y(t) = y(1) * exp(-(t - 1) / RC) for t ≥ 1
    idx2 = t >= 1;  % Indices where t ≥ 1
    y(idx2) = y1 * exp(-(t(idx2) - 1) / RC);

    % --- PLOTTING ---
    figure;
    plot(t, x, 'k--', 'LineWidth', 1.5); hold on; % Input x(t) in dashed black
    plot(t, y, 'b-', 'LineWidth', 2);  % Output y(t) in solid blue

    % Label and title
    xlabel('Time (t)');
    ylabel('Amplitude');
    title(['RC Circuit Response for RC = ', num2str(RC)]);

    % Legend and grid for clarity
    legend('Input x(t)', 'Output y(t)');
    grid on;
end    

