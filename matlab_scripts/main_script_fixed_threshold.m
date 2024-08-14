clear;
close all;
clc;

%% Simulation Parameter Section

%Sampling parameters
fs = 200e3; %Sampling frequency in Hz
duration = 5; %Overall sampling duration in s (duration of the simulation)

%AWGN-Channel parameters
awgn_snr = -15; %SNR after the AWGN-channel in dB

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
time_code_signal = 0.5*(square(-2*pi*1*t, 80) + 1);

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
s_prev = 0;
s_prev2 = 0;
goertzel_segments_results = zeros(1, goertzel_num_segments); %Vector for storing the results of the analysis of the segments 

%Goertzel algorithm execution
for seg = 1:goertzel_num_segments
    %Resetting the Goertzel variables
    s_prev = 0;
    s_prev2 = 0;

    %Calculating the start and end index of the samples to process
    index_start = (seg-1)*goertzel_segment_size + 1;
    index_end = seg*goertzel_segment_size;
    
    %Executing the Goertzel algorithm
    for j = index_start:index_end
        s = dcf77_signal_noise(j) + goertzel_coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    end
    %Calculating the power of the signal in the segment and storing the result in the vector
    goertzel_segments_results(seg) = s_prev2^2 + s_prev^2 - goertzel_coeff * s_prev * s_prev2;
end



%% Detection Section

%Normalize the results
goertzel_segments_results_normalized = goertzel_segments_results / max(goertzel_segments_results);

%Vector for storing the results of the analysis of the segments
dcf77_reconstructed = zeros(1, goertzel_num_segments);

%Looping over the segments and checking if value is below or above the threshold
for seg = 1:goertzel_num_segments
    if goertzel_segments_results_normalized(seg) > 0.15
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
subplot(7,1,1);
plot(t, carrier_signal);
title('Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(7,1,2);
plot(t, time_code_signal);
title('Time-Code signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(7,1,3);
plot(t, dcf77_signal);
title('DCF77 Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(7,1,4);
plot(t, dcf77_signal_noise);
title('DCF77 Signal with noise');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(7,1,5);
stem(t_goertzel_segments_results, goertzel_segments_results);
title('Result Goertzel algorithm');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(7,1,6);
stem(t_goertzel_segments_results, goertzel_segments_results_normalized);
title('Result Goertzel algorithm normalized');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(7,1,7);
stairs(t_goertzel_segments_results, dcf77_reconstructed);
title('Time-Code signal after detector');
xlabel('Time (s)');
ylabel('Amplitude');

%Frequency-Domain plots
figure
subplot(2, 1, 1);
plot(fft_freq, dcf77_signal_fft);
title('Amplitude Spectrum of the DCF77 signal');
xlabel('Frequency (Hz)');
ylabel('20*log_{10}(|DCF77(f)|)');

subplot(2, 1, 2);
plot(fft_freq, dcf77_signal_noise_fft);
title('Amplitude Spectrum of the DCF77 signal after AWGN-Channel');
xlabel('Frequency (Hz)');
ylabel('20*log_{10}(|DCF77-AWGN(f)|)');