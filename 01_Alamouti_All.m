% Matlab MIMO Simulations - 
% Part 1: Alamouti 2*1
% Editted by Ali B. Mar. 2025

% =======================================================
% ======================= Clear All =========================
clc; clear; close;

% =======================================================
% ====================== Parameters =======================
M = 2;                                                   % Number of transmit antennas (Alamouti: 2)
N = 1;                                                   % Number of receive antennas
snr_dB = 0:2:30;                                  % SNR range in dB
snr_lin = 10.^(snr_dB/10);                   % SNR linear scale
N_symbols = 1e6;                               % Number of symbols per modulation

% Modulation settings: [ModOrder, ModType, LegendLabel]
mod_list = {
    2,   'psk',  'BPSK';
    4,   'psk',  'QPSK';
    8,   'psk',  '8PSK';
    16,  'qam',  '16QAM';
    64,  'qam',  '64QAM';
    256, 'qam',  '256QAM'
};
ber = zeros(size(mod_list,1), length(snr_dB));    % BER result storage

% =======================================================
% ======================= Main Loop =======================

for m_idx = 1:size(mod_list,1)
    Mod_order = mod_list{m_idx,1};
    mod_type = mod_list{m_idx,2};
    k = log2(Mod_order);                    % Bits per symbol
    N_bits = N_symbols * k;                 % Total number of bits
    
    for i = 1:length(snr_dB)
        % Generate random bits
        bits = randi([0 1], N_bits, 1);     % Generate bits as column vector
        
        % Modulation
        if strcmp(mod_type, 'psk')
            % Reshape bits into k-bit groups and convert to integers
            bit_groups = reshape(bits, k, [])';    % Reshape to N_symbols × k matrix
            symbols_int = bi2de(bit_groups, 'left-msb');
            symbols = pskmod(symbols_int, Mod_order, 0, 'gray');
            symbols = symbols / sqrt(mean(abs(symbols).^2)); % Normalize
        else
            symbols = qammod(bits, Mod_order, 'InputType', 'bit', ...
                'UnitAveragePower', true);
        end
        
        % Reshape for Alamouti encoding (ensure even number of symbols)
        symbols = reshape(symbols, 2, []);
        
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
        s_hat = s_hat(:).';
        
        % Channel power normalization
        h_power = sum(abs([h1; h2]).^2, 1);
        h_power = repmat(h_power, 2, 1);
        s_hat = s_hat(:) ./ h_power(:);
        
        % Demodulation
        if strcmp(mod_type, 'psk')
            bits_hat = pskdemod(s_hat, Mod_order, 0, 'gray');
            bits_hat = reshape(de2bi(bits_hat, k, 'left-msb')', [], 1);
        else
            bits_hat = qamdemod(s_hat, Mod_order, 'OutputType', 'bit', ...
                'UnitAveragePower', true);
        end
        
        % Compute BER
        ber(m_idx, i) = sum(bits ~= bits_hat) / length(bits);
    end
end

% =======================================================
% ======================= Plot Results ======================

figure('Position', [100 100 800 600]); % Make figure larger
lineStyles = {'b-o', 'r--s', 'g-.^', 'm:d', 'k-v',  'c--*' };

for m_idx = 1:size(mod_list,1)
    semilogy(snr_dB, ber(m_idx,:), lineStyles{m_idx}, 'LineWidth', 2,'MarkerSize', 10, 'MarkerFaceColor', 'w');   
    hold on;
end

% Grid 
grid on;
set(gca, 'GridLineStyle', ':');    
set(gca, 'MinorGridLineStyle', ':');
grid minor;

% Fonts and Labels
set(gca, 'FontSize', 18);          
xlabel('SNR (dB)', 'FontSize', 20);
ylabel('Bit Error Rate (BER)', 'FontSize', 20);
title('Alamouti 2×1 STBC Performance over Rayleigh Fading', 'FontSize', 20, 'FontWeight', 'bold');

% Legend
lgd = legend(mod_list{:,3}, 'Location', 'southwest', 'FontSize', 18);
set(lgd, 'Box', 'on');            

% Set Axis 
ylim([1e-6 1]);
xlim([min(snr_dB) max(snr_dB)]);

% Add Box
box on;