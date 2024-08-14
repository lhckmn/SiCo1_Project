clear;
close all;
clc;

% Parameters
fs = 1000; % Sampling frequency in Hz
duration = 1; % Duration of the signal in seconds
duty_cycle = 0.8; % Duty cycle (80%)

% Time vector
t = 0:1/fs:duration;

% Number of samples
num_samples = length(t);

% Calculate the number of samples for the pulse width and gap
pulse_width_samples = round(duty_cycle * num_samples);
gap_samples = num_samples - pulse_width_samples;

% Create the rectangular pulse with the gap at the beginning
rect_signal = [zeros(1, gap_samples), ones(1, pulse_width_samples)];

% Truncate or extend the signal to match the duration exactly
rect_signal = rect_signal(1:num_samples);

% Plot the signal
figure;
plot(t, rect_signal, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Rectangular Pulse with 80% Duty Cycle and Gap at Beginning');
grid on;

% Optionally, you can save the signal to a file
audiowrite('rect_signal.wav', rect_signal, fs);
