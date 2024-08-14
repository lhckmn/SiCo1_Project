clear;
close all;
clc;


% Parameters
fs = 1000;            % Sampling frequency
f = 100;              % Carrier frequency
duration = 4;         % Duration of the signal in seconds
N = fs * duration;    % Total number of samples
t = (0:N-1)/fs;       % Time vector

% Create a signal with amplitude switching
signal = [sin(2*pi*f*t(1:N/4)), 0.5*sin(2*pi*f*t(N/4+1:N/2)), ...
          sin(2*pi*f*t(N/2+1:3*N/4)), 0.5*sin(2*pi*f*t(3*N/4+1:end))];

% Initialize Goertzel algorithm parameters
k = round(f*duration);            % Frequency bin
omega = 2 * pi * k / N;
coeff = 2 * cos(omega);
s_prev = 0;
s_prev2 = 0;

% Segment size for analysis (e.g., 1 second)
segment_size = fs;

% Initialize array to store power of the carrier frequency
power = zeros(1, duration);

for i = 1:duration
    s_prev = 0;
    s_prev2 = 0;
    for j = 1:segment_size
        s = signal((i-1)*segment_size + j) + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    end
    % Calculate power
    power(i) = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;
end

% Normalize power to see the variations clearly
normalized_power = power / max(power);

% Plot the results
figure;
subplot(2,1,1);
plot(t, signal);
title('Signal with Amplitude Switching');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
stem(1:duration, normalized_power);
title('Normalized Power of Carrier Frequency over Time');
xlabel('Time (s)');
ylabel('Normalized Power');
