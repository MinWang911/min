clear;
clc;
close all;

rng(1);

%% 1. OFDM configuration

cfg.N_fft = 64;
cfg.N_cp = 16;

cfg.N_data = 48;
cfg.N_pilot = 4;

cfg.ModOrder = 4;       % QPSK
cfg.numSymbols = 200;

cfg.Fs = 20e6;

% IEEE 802.11a subcarrier numbers
data_idx_80211 = [ ...
    -26:-22, ...
    -20:-8, ...
    -6:-1, ...
    1:6, ...
    8:20, ...
    22:26];

pilot_idx_80211 = ...
    [-21, -7, 7, 21];

% Convert to MATLAB fftshift indexes
to_matlab_idx = ...
    @(x) x + cfg.N_fft/2 + 1;

cfg.data_idx = ...
    to_matlab_idx(data_idx_80211);

cfg.pilot_idx = ...
    to_matlab_idx(pilot_idx_80211);

%% 2. Run transmitter

[tx_signal, tx_bits] = ...
    wifi_tx_pipeline(cfg);

%% 3. Multipath channel

h = [ ...
    0.9, ...
    0.4*exp(1j*pi/3), ...
    0.2*exp(-1j*pi/4)];

h = h / sqrt(sum(abs(h).^2));

tx_faded = filter( ...
    h, ...
    1, ...
    tx_signal);

%% 4. Add AWGN

SNRdB = 20;

signalPower = ...
    mean(abs(tx_faded).^2);

noisePower = ...
    signalPower / 10^(SNRdB/10);

noise = sqrt(noisePower/2) .* ...
    (randn(size(tx_faded)) + ...
     1j*randn(size(tx_faded)));

rx_signal = ...
    tx_faded + noise;

%% 5. Run receiver

[rx_bits, ber, rx_data_symbols, debug] = ...
    wifi_rx_pipeline( ...
        rx_signal, ...
        cfg, ...
        tx_bits);

%% 6. Display receiver results

numErrors = ...
    sum(rx_bits ~= tx_bits);

accuracy = ...
    100 * (1 - ber);

fprintf('\nOFDM Receiver Result\n');
fprintf('Modulation order = %d\n', cfg.ModOrder);
fprintf('SNR = %.1f dB\n', SNRdB);
fprintf('Total bits = %d\n', length(tx_bits));
fprintf('Recovered bits = %d\n', length(rx_bits));
fprintf('Bit errors = %d\n', numErrors);
fprintf('BER = %.6f\n', ber);
fprintf('Recovery accuracy = %.3f %%\n', accuracy);

%% 7. Receiver front-end result

symbolIndex = 1;

figure;

subplot(3,1,1);

plot(real( ...
    debug.rx_time_no_cp(:, symbolIndex)));

grid on;

xlabel('Sample index');
ylabel('Real amplitude');

title('After Cyclic Prefix Removal');

subplot(3,1,2);

stem( ...
    0:cfg.N_fft-1, ...
    abs(debug.R_unshifted(:, symbolIndex)), ...
    'filled');

grid on;

xlabel('FFT bin index');
ylabel('Magnitude');

title('After FFT');

subplot(3,1,3);

subcarrierAxis = ...
    -cfg.N_fft/2 : cfg.N_fft/2-1;

stem( ...
    subcarrierAxis, ...
    abs(debug.R_shifted(:, symbolIndex)), ...
    'filled');

grid on;

xlabel('Subcarrier index');
ylabel('Magnitude');

title('After FFT Shift');

%% 8. Transmitted and recovered bit comparison

numBitsToShow = ...
    min(100, length(tx_bits));

figure;

stairs( ...
    1:numBitsToShow, ...
    tx_bits(1:numBitsToShow), ...
    'LineWidth', 1.5);

hold on;

stairs( ...
    1:numBitsToShow, ...
    rx_bits(1:numBitsToShow), ...
    '--', ...
    'LineWidth', 1.5);

grid on;

ylim([-0.2 1.2]);

xlabel('Bit index');
ylabel('Bit value');

legend( ...
    'Transmitted bits', ...
    'Recovered bits');

title(sprintf( ...
    'Bit Recovery at %.1f dB, BER = %.4g', ...
    SNRdB, ber));