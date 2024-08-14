clear;
close all;
clc;




% Parameters
fs = 1000;            % Sampling frequency
f = 100;              % Carrier frequency
total_duration = 8;   % Total duration of the signal in seconds
N = fs * total_duration;  % Total number of samples
t = (0:N-1)/fs;       % Time vector

% Create a signal with amplitude switching every second
signal = [sin(2*pi*f*t(1:N/4)), 0.5*sin(2*pi*f*t(N/4+1:N/2)), ...
          sin(2*pi*f*t(N/2+1:3*N/4)), 0.5*sin(2*pi*f*t(3*N/4+1:end))];

% Parameters for Goertzel algorithm
segment_duration = 1;   % Duration of each segment in seconds
segment_size = fs * segment_duration; % Segment size in samples
num_segments = total_duration / segment_duration; % Number of segments

% Initialize Goertzel algorithm parameters
k = round(f*segment_duration);  % Frequency bin for the segment
omega = 2 * pi * k / segment_size;
coeff = 2 * cos(omega);

% Initialize array to store power of the carrier frequency
power = zeros(1, num_segments);

% Loop through each segment and apply Goertzel algorithm
for seg = 1:num_segments
    s_prev = 0;
    s_prev2 = 0;
    start_idx = (seg-1)*segment_size + 1;
    end_idx = seg*segment_size;
    for j = start_idx:end_idx
        s = signal(j) + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    end
    % Calculate power
    power(seg) = s_prev2^2 + s_prev^2 - coeff * s_prev * s_prev2;
end

% Normalize power to see the variations clearly
normalized_power = power / max(power);

% Parameters for adaptive threshold
window_size = 3; % Size of the moving window
alpha = 1.5;    % Multiplier for standard deviation to define the threshold

% Initialize array to store adaptive thresholds
adaptive_threshold = zeros(1, num_segments);

% Calculate adaptive thresholds
for i = window_size:num_segments
    window_data = normalized_power(i-window_size+1:i);
    moving_avg = mean(window_data);
    moving_std = std(window_data);
    adaptive_threshold(i) = moving_avg - alpha * moving_std;
end

% Detect changes using adaptive threshold
changes_adaptive = find(normalized_power < adaptive_threshold);

% Plot the results
figure;
subplot(3,1,1);
plot(t, signal);
title('Signal with Amplitude Switching');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
stem((1:num_segments)*segment_duration, normalized_power);
title('Normalized Power of Carrier Frequency over Time');
xlabel('Time (s)');
ylabel('Normalized Power');

subplot(3,1,3);
stem((1:num_segments)*segment_duration, normalized_power, 'b');
hold on;
plot((1:num_segments)*segment_duration, adaptive_threshold, 'r--');
stem(changes_adaptive * segment_duration, normalized_power(changes_adaptive), 'g');
title('Detected Changes using Adaptive Threshold');
xlabel('Time (s)');
ylabel('Normalized Power');
legend('Power', 'Adaptive Threshold', 'Detected Changes');
hold off;
