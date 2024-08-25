clear;
close all;
clc;

%% Simulation Parameter Section

%Sampling parameters
fs = 24e3; %Sampling frequency in Hz
duration = 30; %Overall sampling duration in s (duration of the simulation)

%AWGN-Channel parameters
awgn_snr = -5; %SNR after the AWGN-channel in dB

%Goertzel macro parameters
goertzel_analyzation_freq = 5.5e3; %Analyzation frequency of the Goertzel algorithm in Hz
goertzel_segment_duration = 0.01; %Duration of the Goertzel algorithm segment in s

%Detector parameters
det_avg_time = 10; %Time of the Exponential Moving Average in s



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

%Concatenating the symbols after each other to create the time code signal
time_code_signal = [];
for i = 1:ceil(duration)
    if mod(i, 2) == 0 %Alternating the symbols
        time_code_signal = [time_code_signal, rect_duty_90];
    else
        time_code_signal = [time_code_signal, rect_duty_80];
    end
end
time_code_signal = time_code_signal(1:N_samples-1); %Cutting to size
time_code_signal = [zeros(1, 1), time_code_signal]; %Inserting a zero at the beginning

%Creating the DCF77 signal by using Amplitude Modulation (carrier is reduced to 15% during the gap)
dcf77_signal = (0.85*time_code_signal + 0.15) .* carrier_signal;

%Add noise by passing the DCF77 signal through an AWGN-Channel
dcf77_signal_noise = awgn(dcf77_signal, awgn_snr, 'measured');

%Verifying SNR by measuring again
snr_measured = snr(dcf77_signal, dcf77_signal_noise - dcf77_signal);



%% Goertzel Algoritm Section

%Batch analysis parameters
goertzel_segment_size = fs * goertzel_segment_duration; %Calculating the number of samples for one segment
goertzel_num_segments = duration / goertzel_segment_duration; %Calculating the number of segments
t_goertzel_segments_results = (1:goertzel_num_segments)*goertzel_segment_duration; %Create a time vector for the results for plotting

%Goertzel parameters
goertzel_k = round(goertzel_analyzation_freq * goertzel_segment_duration);
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
    goertzel_segments_magnitudes(seg) = sqrt(Q_1^2 + Q_2^2 - Q_1 * Q_2 * goertzel_coeff);
end



%% Detection Section

%Preparing vectors for the reconstructed signal and the Exponential Moving Average (EMA)
exp_mov_avg = zeros(1, goertzel_num_segments);
detector_threshold = zeros(1, goertzel_num_segments);
dcf77_reconstructed = zeros(1, goertzel_num_segments);

%Calculating the parameter alpha of the EMA
exp_mov_avg_alpha = 2 / (det_avg_time * (1 / goertzel_segment_duration) + 1);

for seg = 1:goertzel_num_segments
    %Calculating the Exponential Moving Average
    if seg == 1 %Check if it is the first segment
        exp_mov_avg(seg) = exp_mov_avg_alpha * goertzel_segments_magnitudes(seg); %Adapt calulation as there is no element at index -1
    else
        exp_mov_avg(seg) = (1 - exp_mov_avg_alpha) * exp_mov_avg(seg - 1) + exp_mov_avg_alpha * goertzel_segments_magnitudes(seg); %Calculate EMA normally
    end

    %Calculating the current threshold
    detector_threshold(seg) = 0.575 * (1 / 0.8725) * exp_mov_avg(seg);

    %Deciding by comparing sample to threshold
    if goertzel_segments_magnitudes(seg) > detector_threshold(seg)
        dcf77_reconstructed(seg) = 1;
    else
        dcf77_reconstructed(seg) = 0;
    end
end



%% FFT and Plot Section

%Getting the number of samples and preparing a vector for the frequency axis
fft_N = length(t); %Number of samples
fft_freq = (-fft_N/2:fft_N/2-1)*(fs/fft_N); %Frequency vector

%Creating a FFT of the DCF77 signal without noise
dcf77_signal_fft = fft(dcf77_signal); %FFT of DCF77 signal
dcf77_signal_fft = fftshift(dcf77_signal_fft); %Shifting FFT
dcf77_signal_fft = dcf77_signal_fft / fft_N; %Adjusting FFT gain
dcf77_signal_fft = 20*log10(abs(dcf77_signal_fft)); %Converting to logarithmic

%Creating a FFT of the DCF77 signal with noise
dcf77_signal_noise_fft = fft(dcf77_signal_noise); %FFT of DCF77 signal with noise
dcf77_signal_noise_fft = fftshift(dcf77_signal_noise_fft); %Shifting FFT
dcf77_signal_noise_fft = dcf77_signal_noise_fft / fft_N; %Adjusting FFT gain
dcf77_signal_noise_fft = 20*log10(abs(dcf77_signal_noise_fft)); %Converting to logarithmic


%Frequency domain plots
figure
subplot(2, 1, 1);
plot(fft_freq, dcf77_signal_fft);
title('Amplitude Spectrum of the DCF77 Signal');
xlabel('Frequency (Hz)');
ylabel('20*log_{10}(|DCF77(f)|)');

subplot(2, 1, 2);
plot(fft_freq, dcf77_signal_noise_fft);
title('Amplitude Spectrum of the DCF77 Signal after AWGN-Channel');
xlabel('Frequency (Hz)');
ylabel('20*log_{10}(|DCF77_{AWGN}(f)|)');


%Time domain plots
figure

subplot(5,1,1);
plot(t, time_code_signal);
title('Time-Code Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,2);
plot(t, dcf77_signal);
title('DCF77 Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,3);
plot(t, dcf77_signal_noise);
title('DCF77 Signal with noise');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(5,1,4);
stem(t_goertzel_segments_results, goertzel_segments_magnitudes);
hold on;
stairs(t_goertzel_segments_results, detector_threshold, 'Color', 'red');
title('Result Goertzel algorithm and Detector Threshold');
xlabel('Time (s)');
ylabel('Magnitude');
hold off;

% subplot(5,1,5);
% stairs(t_goertzel_segments_results, exp_mov_avg);
% title('EMA of the detector');
% xlabel('Time (s)');
% ylabel('Value');

subplot(5,1,5);
stairs(t_goertzel_segments_results, dcf77_reconstructed);
title('Time-Code Signal after detector');
xlabel('Time (s)');
ylabel('Value');