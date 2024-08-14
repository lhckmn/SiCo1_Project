clear;
close all;
clc;

%% Simulation Parameter Section

%Sampling parameters
fs = 200e3; %Sampling frequency in Hz
duration = 24; %Overall sampling duration in s (duration of the simulation)

%AWGN-Channel parameters
awgn_snr = -10; %SNR after the AWGN-channel in dB

%Goertzel Macro parameters
goertzel_segment_duration = 0.01; %Duration of the Goertzel algorithm segment in s

%Detector parameters
det_threshold = 0.15; %Threshold for the detector

%% Signal Creation Section

%Creating a time vector for the duration of the simulation
t = 0:1/fs:duration;

%Creating a variable with the number of samples
N_samples = length(t);

%Carrier signal (cos at 77.5kHz)
fc = 77.5e3; %Storing the carrier frequency of the DCF77 transmitter in Hz (77.5kHz)
carrier_signal = cos(2*pi*fc*t);

%Time signal (square wave with frequency of 1Hz, 80% duty cycle, and a range from 0 to 1)
one_samples_80 = round(0.8 * fs);
one_samples_90 = round(0.9 * fs);

rect_duty_80 = [zeros(1, fs - one_samples_80), ones(1, one_samples_80)]; %Rect, duty-cycle 80%, starting at 0
rect_duty_90 = [zeros(1, fs - one_samples_90), ones(1, one_samples_90)]; %Rect, duty-cycle 90%, starting at 0

time_code_signal = [];
for i = 1:ceil(duration)
    if mod(i, 3) == 0
        time_code_signal = [time_code_signal, rect_duty_90];
    else
        time_code_signal = [time_code_signal, rect_duty_80];
    end
end
time_code_signal = time_code_signal(1:N_samples-1);
time_code_signal = [zeros(1, 1), time_code_signal];

%Creating the DCF77 signal by using Amplitude Modulation (carrier is reduced to 15% during the gap)
dcf77_signal = (0.85*time_code_signal + 0.15) .* carrier_signal;



%% AWGN-Channel Section

%Add noise by passing the DCF77 signal through an AWGN-Channel
dcf77_signal_noise = awgn(dcf77_signal, awgn_snr, 'measured');



%% Goertzel Section

%Batch analysis parameters
goertzel_segment_size = fs * goertzel_segment_duration; %Calculating the number of samples for one segment
goertzel_num_segments = duration / goertzel_segment_duration; %Calculating the number of segments
t_goertzel_segments_results = (1:goertzel_num_segments)*goertzel_segment_duration; %Create a time vector for the results for plotting

%Goertzel parameters
goertzel_k = round(fc * goertzel_segment_duration);
goertzel_omega = 2 * pi * goertzel_k / goertzel_segment_size;
goertzel_coeff = 2 * cos(goertzel_omega);

%Initializing Goertzel variables and a vector to store the results
Q_0 = 0;
Q_1 = 0;
Q_2 = 0;
goertzel_segments_magnitudes = zeros(1, goertzel_num_segments); %Vector for storing the results of the analysis of the segments 

%Goertzel algorithm segment processing
for seg = 1:goertzel_num_segments
    %Resetting the Goertzel variables
    Q_0 = 0;
    Q_1 = 0;
    Q_2 = 0;

    %Calculating the start and end index of the samples to process
    index_start = (seg-1)*goertzel_segment_size + 1;
    index_end = seg*goertzel_segment_size;
    
    %Executing the Goertzel algorithm
    for j = index_start:index_end
        Q_0 = Q_1 * goertzel_coeff - Q_2 + dcf77_signal_noise(j);
        Q_2 = Q_1;
        Q_1 = Q_0;
    end
    %Calculating the squared magnitude of fc in the segment and storing the result in the vector of segment results
    goertzel_segments_magnitudes(seg) = Q_1^2 + Q_2^2 - Q_1 * Q_2 * goertzel_coeff;
end



%% Detection Section

%Preparing vectors for the reconstructed signal and a rolling average 
rolling_avg = zeros(1, goertzel_num_segments);
goertzel_segments_magnitudes_normalized = zeros(1, goertzel_num_segments);
dcf77_reconstructed = zeros(1, goertzel_num_segments);

for seg = 1:goertzel_num_segments
    %Calculating the rolling average
    if seg == 1
        rolling_avg = goertzel_segments_magnitudes(seg); %Initialize with current value on first segment
    else
        rolling_avg(seg) = 0.999 * rolling_avg(seg - 1) + 0.001 * goertzel_segments_magnitudes(seg); %Calculate a weighted average otherwise
    end

    %Normalize current value with average
    goertzel_segments_magnitudes_normalized(seg) = goertzel_segments_magnitudes(seg) / rolling_avg(seg);

    if goertzel_segments_magnitudes_normalized(seg) > 0.2
        dcf77_reconstructed(seg) = 1;
    else
        dcf77_reconstructed(seg) = 0;
    end
end



%% FFT Section

fft_N = length(t); %Number of samples
fft_freq = (-fft_N/2:fft_N/2-1)*(fs/fft_N); %Frequency vector


dcf77_signal_fft = fft(dcf77_signal); %FFT of DCF77 signal
dcf77_signal_fft = fftshift(dcf77_signal_fft); %Shifting FFT
dcf77_signal_fft = dcf77_signal_fft / fft_N; %Adjusting FFT gain
dcf77_signal_fft = 20*log10(abs(dcf77_signal_fft)); %Converting to logarithmic

dcf77_signal_noise_fft = fft(dcf77_signal_noise); %FFT of DCF77 signal with noise
dcf77_signal_noise_fft = fftshift(dcf77_signal_noise_fft); %Shifting FFT
dcf77_signal_noise_fft = dcf77_signal_noise_fft / fft_N; %Adjusting FFT gain
dcf77_signal_noise_fft = 20*log10(abs(dcf77_signal_noise_fft)); %Converting to logarithmic


%% Plot Section

%Time-Domain plots
figure
subplot(8,1,1);
plot(t, carrier_signal);
title('Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(8,1,2);
plot(t, time_code_signal);
title('Time-Code signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(8,1,3);
plot(t, dcf77_signal);
title('DCF77 Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(8,1,4);
plot(t, dcf77_signal_noise);
title('DCF77 Signal with noise');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(8,1,5);
stem(t_goertzel_segments_results, goertzel_segments_magnitudes);
title('Result Goertzel algorithm');
xlabel('Time (s)');
ylabel('Magnitude^2');

subplot(8,1,6);
stairs(t_goertzel_segments_results, rolling_avg);
title('Rolling average of the detector');
xlabel('Time (s)');
ylabel('Value');

subplot(8,1,7);
stem(t_goertzel_segments_results, goertzel_segments_magnitudes_normalized);
title('Result Goertzel algorithm normalized');
xlabel('Time (s)');
ylabel('Value');

subplot(8,1,8);
stairs(t_goertzel_segments_results, dcf77_reconstructed);
title('Time-Code signal after detector');
xlabel('Time (s)');
ylabel('Value');


% %Frequency-Domain plots
% figure
% subplot(2, 1, 1);
% plot(fft_freq, dcf77_signal_fft);
% title('Amplitude Spectrum of the DCF77 signal');
% xlabel('Frequency (Hz)');
% ylabel('20*log_{10}(|DCF77(f)|)');
% 
% subplot(2, 1, 2);
% plot(fft_freq, dcf77_signal_noise_fft);
% title('Amplitude Spectrum of the DCF77 signal after AWGN-Channel');
% xlabel('Frequency (Hz)');
% ylabel('20*log_{10}(|DCF77-AWGN(f)|)');