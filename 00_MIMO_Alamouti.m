% Matlab MIMO Simulations - 
% Part 1: Alamouti 2*1
% Editted by Ali B. Mar. 2025

clc; clear; close;

% ====================================================
% ===================== Parameters =====================

M = 2;                                       % Number of transmit antennas 
N = 1;                                       % Number of receive antennas
snr_dB = 1:1:20;                      % SNR range in dB
snr_lin = 10.^(snr_dB/10);        % SNR linear scale
N_bits = 1e6;                           % Total number of bits
ber = zeros(size(snr_dB));      % BER result storage               

% ====================================================
% ===================== Main Loop ======================

for i=1:length(snr_dB)
    
    bits = randi([0 1], 1, N_bits);          % Generate random bits
    symbols = pskmod(bits, 2, 0);       % BPSK with 0 is the phase offset
    symbols = reshape(symbols, M, []);

    % Rayleigh flat fading channels from each TX antenna
    h1 = (randn(1, size(symbols,2)) + 1j*randn(1, size(symbols,2))) / sqrt(2);
    h2 = (randn(1, size(symbols,2)) + 1j*randn(1, size(symbols,2))) / sqrt(2);

    % Alamouti Encoding
    x1 = symbols(1,:);
    x2 = symbols(2,:);

    tx1 = [x1; -conj(x2)];   % Antenna 1
    tx2 = [x2;  conj(x1)];   % Antenna 2

    % AWGN noise
    noise_var = 1 / snr_lin(i);
    n1 = sqrt(noise_var/2) * (randn(1, size(symbols,2)) + 1j*randn(1, size(symbols,2)));
    n2 = sqrt(noise_var/2) * (randn(1, size(symbols,2)) + 1j*randn(1, size(symbols,2)));

    % Received signals
    r1 = h1 .* tx1(1,:) + h2 .* tx2(1,:) + n1;
    r2 = h1 .* tx1(2,:) + h2 .* tx2(2,:) + n2;
    % Alamouti Decoding
    s1_hat = conj(h1).*r1 + h2.*conj(r2);
    s2_hat = conj(h2).*r1 - h1.*conj(r2);

    % Stack decoded symbols
    s_hat = [s1_hat; s2_hat];

    % Demodulate using pskdemod
    bits_hat = pskdemod(s_hat, 2, 0);
    bits_rx = reshape(bits_hat, 1, []);

    % Compute BER
    ber(i) = sum(bits_rx ~= bits(1:2*size(symbols,2))) / (2 * size(symbols,2));

end

% =======================================================
% ======================== Plot results ======================
figure;
semilogy(snr_dB, ber, 'b-o', 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 18);          
xlabel('SNR (dB)', 'FontSize', 20);
ylabel('Bit Error Rate (BER)', 'FontSize', 20);
title('Alamouti 2×1 STBC Performance for BPSK over Rayleigh Fading', 'FontSize', 20, 'FontWeight', 'bold');
legend('M = 2, N = 1', 'FontSize', 22);