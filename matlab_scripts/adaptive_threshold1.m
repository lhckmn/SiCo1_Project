clear;
close all;
clc;



% Generate a longer sample signal
rng(0);  % For reproducibility
N = 1000;  % Length of the signal
signal = [randn(1, 300) * 0.1, randn(1, 200) * 0.1 + 1, randn(1, 300) * 0.1, randn(1, 200) * 0.1 + 1];

% Initialize parameters
threshold = 0.5;  % Initial threshold
output = zeros(1, length(signal));  % Initialize output array
alpha = 0.001;  % Learning rate for adaptive threshold

% Adaptive thresholding
for i = 1:length(signal)
    % Determine output based on current threshold
    if signal(i) > threshold
        output(i) = 1;
    else
        output(i) = 0;
    end
    
    % Update the threshold based on the current sample
    % If output is 1, threshold should move closer to the current sample if it is 1
    % If output is 0, threshold should move closer to the current sample if it is 0
    if output(i) == 1
        threshold = threshold + alpha * (signal(i) - threshold);
    else
        threshold = threshold - alpha * (threshold - signal(i));
    end
end

% Display results
disp('Final Threshold:');
disp(threshold);

% Plot the signal, threshold, and output
figure;
subplot(3,1,1);
plot(signal, '-o');
title('Input Signal');
ylabel('Amplitude');
xlabel('Sample');

subplot(3,1,2);
plot(output, '-o');
title('Output Signal');
ylabel('Binary Output');
xlabel('Sample');

subplot(3,1,3);
plot(1:length(signal), threshold * ones(1,length(signal)), '--');
title('Adaptive Threshold');
ylabel('Threshold');
xlabel('Sample');
