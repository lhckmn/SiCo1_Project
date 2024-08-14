clear;
close all;
clc;


% Parameters
fs = 8000;  % Sampling frequency
f = 1000;   % Frequency to detect
N = 256;    % Number of samples
n = 0:N-1;  % Sample index

% Generate a test signal (sine wave at 1 kHz)
x = sin(2 * pi * f * n / fs);

% Apply the Goertzel Algorithm
k = round((f / fs) * N);  % Bin index for target frequency
result = goertzel(x, k);

% Display result
disp(['Magnitude at ', num2str(f), ' Hz: ', num2str(abs(result))]);

function magnitude = goertzel(x, k)
    N = length(x);
    omega = 2 * pi * k / N;
    cos_omega = cos(omega);
    sin_omega = sin(omega);
    coeff = 2 * cos_omega;

    Q1 = 0;
    Q2 = 0;
    for n = 1:N
        Q0 = x(n) + coeff * Q1 - Q2;
        Q2 = Q1;
        Q1 = Q0;
    end

    real_part = (Q1 - Q2 * cos_omega);
    imag_part = (Q2 * sin_omega);

    magnitude = real_part + 1i * imag_part;
end
