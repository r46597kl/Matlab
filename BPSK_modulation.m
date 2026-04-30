% BPSK Modulation
% This script implements Binary Phase Shift Keying (BPSK) modulation

clear all;
close all;
clc;

% Parameters
Eb = 1;                          % Energy per bit
N0 = 0.1;                        % Noise power spectral density
fc = 1e9;                        % Carrier frequency (1 GHz)
fs = 10 * fc;                    % Sampling frequency
Ts = 1/fs;                       % Sampling period
Tb = 1e-8;                       % Bit duration
Nsamp = fs * Tb;                 % Samples per bit

% Generate random bits
Nbits = 100;
bits = randi([0 1], 1, Nbits);

% Map bits to symbols (+1 for bit 1, -1 for bit 0)
symbols = 2 * bits - 1;

% Upsample the symbols
upsampled = kron(symbols, ones(1, Nsamp));

% Generate carrier signal
t = (0:length(upsampled)-1) * Ts;
carrier = sqrt(2 * Eb / Tb) * cos(2 * pi * fc * t);

% Modulate signal
modulated = upsampled .* carrier;

% Add AWGN noise
noise = sqrt(N0/2) * randn(size(modulated));
received = modulated + noise;

% Demodulate signal
demodulated = received .* carrier;

% Integrate and dump
integrated = zeros(1, Nbits);
for i = 1:Nbits
    idx = (i-1)*Nsamp + 1:i*Nsamp;
    integrated(i) = sum(demodulated(idx)) / Nsamp;
end

% Threshold detection
threshold = 0;
detected_symbols = sign(integrated);
detected_bits = (detected_symbols + 1) / 2;

% Calculate BER
errors = sum(abs(bits - detected_bits));
BER = errors / Nbits;

% Display results
fprintf('BPSK Modulation Results\n');
fprintf('=======================\n');
fprintf('Number of bits: %d\n', Nbits);
fprintf('Errors: %d\n', errors);
fprintf('Bit Error Rate: %f\n', BER);

% Plot results
figure('Position', [100 100 1200 400]);

subplot(1, 3, 1);
plot(t(1:1000), modulated(1:1000));
xlabel('Time (s)');
ylabel('Amplitude');
title('BPSK Modulated Signal');
grid on;

subplot(1, 3, 2);
plot(t(1:1000), received(1:1000));
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal (with noise)');
grid on;

subplot(1, 3, 3);
stem(1:Nbits, bits, 'b', 'filled');
hold on;
stem(1:Nbits, detected_bits, 'r--');
xlabel('Bit index');
ylabel('Bit value');
title('Original vs Detected Bits');
legend('Original', 'Detected');
grid on;
