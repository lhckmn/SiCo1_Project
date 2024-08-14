clear;
close all;
clc;



% Define the time vector
t = 0:0.01:10; % from 0 to 10 seconds with a step of 0.01 seconds

% Define the ramp signal
ramp1 = t(t<=2);         % Ramp from 0 to 2 seconds
ramp2 = 2*ones(size(t(t>2 & t<=4))); % Constant from 2 to 4 seconds
ramp3 = (t(t>4 & t<=6)-4) + 2; % Ramp from 4 to 6 seconds starting from 2
ramp4 = 4*ones(size(t(t>6 & t<=8))); % Constant from 6 to 8 seconds
ramp5 = (t(t>8 & t<=10)-8) + 4; % Ramp from 8 to 10 seconds starting from 4

% Concatenate all parts to form the complete signal
signal = [ramp1, ramp2, ramp3, ramp4, ramp5];

% Calculate the derivative using numerical differentiation
dt = t(2) - t(1); % time step
derivative = diff(signal) / dt; % numerical derivative
t_derivative = t(1:end-1); % time vector for the derivative

% Plot the original signal and its derivative
figure;

% Plot original signal
subplot(2,1,1);
plot(t, signal, 'b', 'LineWidth', 1.5);
title('Original Signal (Series of Ramps)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot derivative of the signal
subplot(2,1,2);
plot(t_derivative, derivative, 'r', 'LineWidth', 1.5);
title('Derivative of the Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Adjust plot spacing
tight_layout();
