% ELCT 222 --- Computer Assignment 2 
% Problem 2(b): MATLAB convolution of two rectangular pulses

% Define time range 
t = -1:0.01:10;

% Define the four test cases: [A, B, T1, T2]
cases = [
    1, 1, 1, 1;
    1, 1, 2, 2;
    -1, -1, 1, 2;
    1, -1, 2, 1;
];

% Loop over all 4 cases
for i = 1:4
    % Extract parameters
    A = cases(i, 1);
    B = cases(i, 2);
    T1 = cases(i, 3);
    T2 = cases(i, 4);

    % Define x(t) = A * [u(t) - u(t - T1)]
    x = A * ((t >= 0) & (t < T1));

    % Define h(t) = B * [u(t) - u(t - T2)]
    h = B * ((t >= 0) & (t < T2));

    % Perform convolution
    y = conv(x, h) * 0.01;  % Multiply by dt to approximate integral

    % Time vector for convolution result
    t_y = 2 * t(1) + 0.01 * (0:length(y)-1);

    % Plot result
    figure;
    plot(t_y, y, 'b', 'LineWidth', 2);
    xlabel('Time (t)');
    ylabel('y(t) = x(t) * h(t)');
    title(sprintf('Convolution Result: A = %d, B = %d, T_1 = %d, T_2 = %d', A, B, T1, T2));
    grid on;
end