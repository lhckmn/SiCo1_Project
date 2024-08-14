clear;
close all;
clc;

%% Parameter Section

%Sampling parameters
fs = 200e3; %Sampling frequency in Hz
duration = 5; %Sampling duration in s

%AWGN-Channel parameters
awgn_snr = 0; %SNR of the AWGN-channel in dB



%% Signal Creation Section

%Create time vector
t = 0:1/fs:duration;

%Create a variable with the number of samples
N_samples = length(t);

%Carrier signal (cos at 77.5kHz)
fc = 77.5e3; %Carrier frequency of the DCF77 transmitter in Hz (77.5kHz)
carrier_signal = cos(2*pi*fc*t);

%Time signal (square wave with frequency of 1Hz, 80% duty cycle, and a range from 0 to 1)
time_signal = 0.5*(square(-2*pi*1*t, 80) + 1);

%Create the DCF77 signal by using Amplitude Modulation (carrier is reduced to 15% during the gap)
dcf77_signal = (0.85*time_signal + 0.15) .* carrier_signal;



%% AWGN-Channel Section

%Add noise by passing the DCF77 signal through an AWGN-Channel
dcf77_signal_noise = awgn(dcf77_signal, awgn_snr, 'measured');



%% Goertzel Section

%Goertzel parameters
goertzel_N = N_samples;                             % Length of signal
goertzel_k = round(fc * goertzel_N / fs);           % Bin index corresponding to carrier frequency
goertzel_omega = (2*pi*goertzel_k) / goertzel_N;    % Frequency bin
goertzel_coeff = 2*cos(goertzel_omega);             % Goertzel coefficient

%Initialize Goertzel variables
s_prev = 0;
s_prev2 = 0;
goertzel_magnitude_vec = zeros(1, goertzel_N);

%Goertzel algorithm execution
for n = 1:goertzel_N
    s = dcf77_signal_noise(n) + goertzel_coeff * s_prev - s_prev2;
    s_prev2 = s_prev;
    s_prev = s;
    goertzel_magnitude_vec(n) = sqrt(s_prev^2 + s_prev2^2 - goertzel_coeff * s_prev * s_prev2);
end



%% Derivative Section

dt = t(2) - t(1);
derivative = diff(goertzel_magnitude_vec) / dt;
derivative_time = t(1:end-1);



%% FFT Section

fft_N = length(t); %Number of samples
fft_freq = (-fft_N/2:fft_N/2-1)*(fs/fft_N); %Frequency vector


dcf77_signal_fft = fft(dcf77_signal); %FFT of DCF77 signal
dcf77_signal_fft = dcf77_signal_fft / fft_N; %Adjust FFT gain of DCF77 signal
dcf77_signal_fft = fftshift(dcf77_signal_fft); %Shift FFT of DCF77 signal

dcf77_signal_noise_fft = fft(dcf77_signal_noise); %FFT of DCF77 signal
dcf77_signal_noise_fft = dcf77_signal_noise_fft / fft_N;
dcf77_signal_noise_fft = fftshift(dcf77_signal_noise_fft); %Shifted FFT of DCF77 signal



%% Plot Section

%Time-Domain plots
figure
subplot(6,1,1);
plot(t, carrier_signal);
title('Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6,1,2);
plot(t, time_signal);
title('Time-Code Signal (1 Hz Square Wave)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6,1,3);
plot(t, dcf77_signal);
title('DCF77 Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6,1,4);
plot(t, dcf77_signal_noise);
title('DCF77 Signal with noise');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6,1,5);
plot(t, goertzel_magnitude_vec);
title('Result Goertzel algorithm');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(6,1,6);
plot(derivative_time, abs(derivative));
title('Derivative Result Goertzel algorithm');
xlabel('Time (s)');
ylabel('Amplitude');

% %Frequency-Domain plots
% figure
% subplot(2, 1, 1);
% plot(fft_freq, abs(dcf77_signal_fft));
% title('Amplitude Spectrum of the DCF77 signal');
% xlabel('Frequency (Hz)');
% ylabel('|DCF77(f)|');
% 
% subplot(2, 1, 2);
% plot(fft_freq, abs(dcf77_signal_noise_fft));
% title('Amplitude Spectrum of the DCF77 signal after AWGN-Channel');
% xlabel('Frequency (Hz)');
% ylabel('|DCF77-AWGN(f)|');