clc;
clear;
close all;

%% ==================== Parameters ====================
numCarr = 256;
cycPrefLen = 32;
nullIdx = [1:32, 129, 225:256]';
numDataCarr = numCarr - length(nullIdx);
ofdmSymbols = 2000;
fc = 2.4e9;
c = 3e8;
h_uav = 100;
Ptx_dBm = 20;
B = numCarr * 15e3;
NF_dB = 5;
P_N = -174 + 10*log10(B) + NF_dB;
K_rician = 10;
BER_target = 1e-2;
a = 9.61; b = 0.16;
L_los = 1.0;
L_nlos = 20;

%% Modulation Schemes
modSchemes = struct( ...
    'name',  {'BPSK', 'QPSK', '16-QAM', '64-QAM'}, ...
    'order', {2,      4,      16,       64}, ...
    'k',     {1,      2,      4,        6}, ...
    'type',  {'psk',  'psk',  'qam',    'qam'} ...
);
nMod = length(modSchemes);

fprintf('====== BER vs Es/N0 Monte Carlo (Mixed Fading) ======\n');
EsN0_dB_vec = 0:2:40;
numSymPerTrial = 2000;
maxTrials = 100;
R_ref = 500;
P_los_ref = calc_Plos(R_ref, h_uav, a, b);
fprintf('Reference P_LOS = %.4f (at R=%dm)\n', P_los_ref, R_ref);
ber_esn0 = zeros(nMod, length(EsN0_dB_vec));
for s = 1:nMod
    M       = modSchemes(s).order;
    bps     = modSchemes(s).k;
    modType = modSchemes(s).type;
    fprintf('\n--- %s ---\n', modSchemes(s).name);
    for i = 1:length(EsN0_dB_vec)
        EsN0_dB  = EsN0_dB_vec(i);
        EsN0_lin = 10^(EsN0_dB / 10);
        totalErrors = 0;
        totalBits   = 0;
        numBits_trial = numSymPerTrial * bps;
        for trial = 1:maxTrials
            srcBits = randi([0,1], numBits_trial, 1);
            if strcmp(modType, 'psk')
                modOut = pskmod(srcBits, M, 'InputType', 'bit');
            else
                modOut = qammod(srcBits, M, 'InputType', 'bit', 'UnitAveragePower', true);
            end
            if rand < P_los_ref
                h = sqrt(K_rician/(K_rician+1)) + ...
                    sqrt(1/(K_rician+1)/2) * (randn(numSymPerTrial,1) + 1j*randn(numSymPerTrial,1));
            else
                h = (randn(numSymPerTrial,1) + 1j*randn(numSymPerTrial,1)) / sqrt(2);
            end
            txFaded  = h .* modOut;
            noisePow = mean(abs(txFaded).^2) / EsN0_lin;
            noise    = sqrt(noisePow/2) * (randn(size(txFaded)) + 1j*randn(size(txFaded)));
            rxSig    = txFaded + noise;
            rxEq     = rxSig ./ h;
            if strcmp(modType, 'psk')
                rxBits = pskdemod(rxEq, M, 'OutputType', 'bit');
            else
                rxBits = qamdemod(rxEq, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            end
            totalErrors = totalErrors + nnz(srcBits ~= rxBits);
            totalBits   = totalBits + numBits_trial;
        end
        ber_esn0(s, i) = max(totalErrors / totalBits, 1e-8);
        fprintf('Es/N0=%2ddB  BER=%.2e\n', EsN0_dB, ber_esn0(s, i));
    end
end

%% AM Thresholds
am_thresh = zeros(1, nMod);
for s = 1:nMod
    ber_s = ber_esn0(s, :);
    valid = ber_s > 0 & ber_s < 0.5;
    [x_u, ia] = unique(log10(ber_s(valid)), 'stable');
    y_u = EsN0_dB_vec(valid); y_u = y_u(ia);
    am_thresh(s) = interp1(x_u, y_u, log10(BER_target), 'pchip');
end
fprintf('\nAM Thresholds:\n');
for s = 1:nMod
    fprintf('  %-7s : %.2f dB\n', modSchemes(s).name, am_thresh(s));
end

EsN0_fine = 0:0.2:40;
ber_esn0_fine = zeros(nMod, length(EsN0_fine));
for s = 1:nMod
    ber_esn0_fine(s, :) = 10.^interp1(EsN0_dB_vec, log10(ber_esn0(s,:)), EsN0_fine, 'pchip');
end
ber_am_snr = NaN(1, length(EsN0_fine));
am_sel_snr = zeros(1, length(EsN0_fine));
for i = 1:length(EsN0_fine)
    selected = 0;
    for s = nMod:-1:1   % 64-QAM
        if ber_esn0_fine(s, i) < BER_target
            selected = s;
            break;
        end
    end
    if selected > 0
        ber_am_snr(i) = ber_esn0_fine(selected, i);
        am_sel_snr(i) = selected;
    end
end

th_colors = {[1 0.6 0], [0 0 1], [0 0.7 0], [1 0 0]};
th_labels = {'th_0', 'th_1', 'th_2', 'th_3'};
figure('Position', [60 60 1000 650]);
hold on;
mod_colors_snr = {[0 0 1], [0 0.7 0], [1 0 0], [0.8 0 0.8]};
mod_styles_snr = {'-', '-', '-', '-'};
for s = 1:nMod
    plot(EsN0_fine, ber_esn0_fine(s,:), mod_styles_snr{s}, ...
        'Color', mod_colors_snr{s}, 'LineWidth', 1.5, ...
        'DisplayName', modSchemes(s).name);
end
valid_am = ~isnan(ber_am_snr);
plot(EsN0_fine(valid_am), ber_am_snr(valid_am), '-k', 'LineWidth', 2.5, ...
    'DisplayName', 'Adaptive Modulation');
yline(BER_target, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(34.5, BER_target*1.3, 'BER = 10^{-2}', ...
    'FontSize', 10, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');

for s = 1:nMod
    xline(am_thresh(s), ':', 'Color', th_colors{s}, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    text(am_thresh(s), 0.2, sprintf('%s=%.2f dB', th_labels{s}, am_thresh(s)), ...
        'FontSize', 9, 'Color', th_colors{s}, 'Rotation', 90, ...
        'VerticalAlignment', 'bottom');
end
hold off;
set(gca, 'YScale', 'log');
grid on;
xlabel('Es/N0 (dB)', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title('BER vs SNR with AM (Mixed Rician/Rayleigh Fading)', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
xlim([0 35]);
ylim([1e-4 1]);
set(gca, 'FontSize', 11);


bps_vec = [1, 2, 4, 6];
EsN0_lin_fine = 10.^(EsN0_fine / 10);

shannon_cap = log2(1 + EsN0_lin_fine);

tp_effective = zeros(nMod, length(EsN0_fine));
tp_raw       = zeros(nMod, length(EsN0_fine));
for s = 1:nMod
    bps = bps_vec(s);
    for i = 1:length(EsN0_fine)
        tp_raw(s, i) = bps * (1 - ber_esn0_fine(s, i));
        if ber_esn0_fine(s, i) < BER_target
            tp_effective(s, i) = bps;
        end
    end
end

% AM throughput
tp_am = zeros(1, length(EsN0_fine));
for i = 1:length(EsN0_fine)
    if am_sel_snr(i) > 0
        tp_am(i) = bps_vec(am_sel_snr(i));
    end
end
figure('Position', [100 40 1100 700]);
hold on;

if am_thresh(1) > 0
    fill([0 am_thresh(1) am_thresh(1) 0], [0 0 13 13], ...
        [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(am_thresh(1)/2, 5, 'Outage', 'FontSize', 14, 'FontWeight', 'bold', ...
        'Color', [0.6 0 0], 'HorizontalAlignment', 'center');
end

% raw throughput k*(1-BER)
dash_colors = {[0 0 1], [0 0.7 0], [1 0 0], [0.8 0 0.8]};
for s = 1:nMod
    plot(EsN0_fine, tp_raw(s,:), '--', 'Color', dash_colors{s}, ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end

for s = 1:nMod
    plot(EsN0_fine, tp_effective(s,:), '-', 'Color', dash_colors{s}, ...
        'LineWidth', 1.5, 'DisplayName', modSchemes(s).name);
end

plot(EsN0_fine, shannon_cap, '-', 'Color', [0.8 0 0], 'LineWidth', 3, ...
    'DisplayName', 'Shannon Capacity');

plot(EsN0_fine, tp_am, '-k', 'LineWidth', 2.5, ...
    'DisplayName', 'Adaptive Modulation');
for s = 1:nMod
    xline(am_thresh(s), ':', 'Color', th_colors{s}, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    text(am_thresh(s), 11, sprintf('%s=%.2f dB', th_labels{s}, am_thresh(s)), ...
        'FontSize', 9, 'Color', th_colors{s}, 'Rotation', 90, ...
        'VerticalAlignment', 'bottom');
end
hold off;
grid on;
xlabel('Es/N_0 (dB)', 'FontSize', 13);
ylabel('Throughput (bit/s/Hz)', 'FontSize', 13);
title('Throughput & Shannon Capacity vs SNR (Mixed Fading)', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
xlim([0 40]);
ylim([0 13]);
set(gca, 'FontSize', 11);
fprintf('\n====== BER vs Distance (Urban, h_UAV=%dm) ======\n', h_uav);
R_mc  = 0:50:1000;
nDist = length(R_mc);
numMC_los  = 50;
numMC_nlos = 50;
BER_mod = zeros(nMod, nDist);
for ri = 1:nDist
    R = R_mc(ri);
    P_los       = calc_Plos(R, h_uav, a, b);
    PL_los_val  = calc_PL_los(R, h_uav, fc, c, L_los);
    PL_nlos_val = calc_PL_nlos(R, h_uav, fc, c, L_nlos);
    SNR_dB_los  = Ptx_dBm - PL_los_val  - P_N;
    SNR_dB_nlos = Ptx_dBm - PL_nlos_val - P_N;
    for s = 1:nMod
        M       = modSchemes(s).order;
        bps     = modSchemes(s).k;
        modType = modSchemes(s).type;
        nBits   = numDataCarr * bps * ofdmSymbols;

        % LOS MC
        errors_los = 0; bits_los = 0;
        for mc = 1:numMC_los
            h_ch = sqrt(K_rician/(K_rician+1)) + ...
                sqrt(1/(K_rician+1)/2) * ...
                (randn(numDataCarr, ofdmSymbols) + 1j*randn(numDataCarr, ofdmSymbols));
            srcBits = randi([0,1], nBits, 1);
            if strcmp(modType, 'psk')
                modSym = pskmod(srcBits, M, 'InputType', 'bit');
            else
                modSym = qammod(srcBits, M, 'InputType', 'bit', 'UnitAveragePower', true);
            end
            modSym = reshape(modSym, numDataCarr, ofdmSymbols);
            txFaded = ofdmmod(h_ch .* modSym, numCarr, cycPrefLen, nullIdx);
            rxSig   = awgn(txFaded, SNR_dB_los, 'measured');
            rxSym   = ofdmdemod(rxSig, numCarr, cycPrefLen, 0, nullIdx) ./ h_ch;
            if strcmp(modType, 'psk')
                rxBits = pskdemod(rxSym, M, 'OutputType', 'bit');
            else
                rxBits = qamdemod(rxSym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            end
            errors_los = errors_los + nnz(srcBits ~= rxBits(:));
            bits_los   = bits_los + nBits;
        end
        ber_los = errors_los / bits_los;

        % NLOS
        errors_nlos = 0; bits_nlos = 0;
        for mc = 1:numMC_nlos
            h_ch = (randn(numDataCarr, ofdmSymbols) + ...
                1j*randn(numDataCarr, ofdmSymbols)) / sqrt(2);
            srcBits = randi([0,1], nBits, 1);
            if strcmp(modType, 'psk')
                modSym = pskmod(srcBits, M, 'InputType', 'bit');
            else
                modSym = qammod(srcBits, M, 'InputType', 'bit', 'UnitAveragePower', true);
            end
            modSym = reshape(modSym, numDataCarr, ofdmSymbols);
            txFaded = ofdmmod(h_ch .* modSym, numCarr, cycPrefLen, nullIdx);
            rxSig   = awgn(txFaded, SNR_dB_nlos, 'measured');
            rxSym   = ofdmdemod(rxSig, numCarr, cycPrefLen, 0, nullIdx) ./ h_ch;
            if strcmp(modType, 'psk')
                rxBits = pskdemod(rxSym, M, 'OutputType', 'bit');
            else
                rxBits = qamdemod(rxSym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            end
            errors_nlos = errors_nlos + nnz(srcBits ~= rxBits(:));
            bits_nlos   = bits_nlos + nBits;
        end
        ber_nlos = errors_nlos / bits_nlos;
        BER_mod(s, ri) = max(P_los * ber_los + (1 - P_los) * ber_nlos, 1e-7);
    end
    fprintf('R=%4dm  P_LOS=%.3f  | BPSK=%.2e  QPSK=%.2e  16QAM=%.2e  64QAM=%.2e\n', ...
        R, P_los, BER_mod(1,ri), BER_mod(2,ri), BER_mod(3,ri), BER_mod(4,ri));
end

R_fine = 0:1:R_mc(end);
nFine  = length(R_fine);
BER_smooth = zeros(nMod, nFine);
cross_dist = zeros(1, nMod);
for s = 1:nMod
    ber_raw = BER_mod(s, :);
    for i = 2:nDist
        ber_raw(i) = max(ber_raw(i), ber_raw(i-1));
    end
    BER_smooth(s, :) = 10.^interp1(R_mc, log10(ber_raw), R_fine, 'pchip');
end
for s = 1:nMod
    idx = find(BER_smooth(s, :) >= BER_target, 1, 'first');
    if isempty(idx)
        cross_dist(s) = R_fine(end) + 1;
    else
        cross_dist(s) = R_fine(idx);
    end
end
fprintf('\nCrossing distances: 64QAM=%dm  16QAM=%dm  QPSK=%dm  BPSK=%dm\n', ...
    cross_dist(4), cross_dist(3), cross_dist(2), cross_dist(1));

BER_am = NaN(1, nFine);
am_sel = zeros(1, nFine);
for i = 1:nFine
    R = R_fine(i);
    selected = 0;
    for s = nMod:-1:1
        if R < cross_dist(s)
            selected = s;
            break;
        end
    end
    if selected > 0
        BER_am(i) = BER_smooth(selected, i);
        am_sel(i) = selected;
    end
end

mod_colors = {[0 0 1], [0 0.7 0], [1 0 0], [0.8 0 0.8]};
mod_styles = {'-', '--', '-.', ':'};
mod_lw     = [1.5, 1.5, 1.5, 1.5];
ylims      = [1e-7 1];


figure('Position', [80 80 1000 650]);
hold on;

for s = 1:nMod
    plot(R_fine, BER_smooth(s, :), mod_styles{s}, ...
        'Color', mod_colors{s}, 'LineWidth', mod_lw(s), ...
        'DisplayName', modSchemes(s).name);
end
yline(BER_target, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(R_fine(end)-5, BER_target*1.4, 'BER_{target} = 10^{-2}', ...
    'FontSize', 10, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
hold off;
set(gca, 'YScale', 'log');
grid on;
xlabel('Horizontal Distance R (m)', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs Distance (Urban, h_{UAV}=%dm, A2G Fading, Monte Carlo)', h_uav), ...
    'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
xlim([0 R_fine(end)]);
ylim(ylims);
set(gca, 'FontSize', 11);

figure('Position', [120 60 1000 650]);
hold on;

for s = 1:nMod
    plot(R_fine, BER_smooth(s, :), mod_styles{s}, ...
        'Color', mod_colors{s}, 'LineWidth', mod_lw(s), ...
        'DisplayName', modSchemes(s).name);
end

valid = ~isnan(BER_am);
plot(R_fine(valid), BER_am(valid), '-k', 'LineWidth', 2.5, ...
    'DisplayName', 'Adaptive Modulation');

yline(BER_target, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(R_fine(end)-5, BER_target*1.4, 'BER_{target} = 10^{-2}', ...
    'FontSize', 10, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');

cross_names = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};
for s = [4, 3, 2, 1]
    cd = cross_dist(s);
    if cd <= R_fine(end)
        xline(cd, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'DisplayName', sprintf('%s=%dm', cross_names{s}, round(cd)));
    end
end

top_y = ylims(2) * 0.6;
cd_sorted = sort(cross_dist);
edges = [0, cd_sorted, R_fine(end)];
[~, sort_idx] = sort(cross_dist);
label_order = cell(1, nMod+1);
rev_names = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};
for k = 1:nMod
    label_order{k} = rev_names{sort_idx(k)};
end
label_order{nMod+1} = 'Outage';
valid_edges = edges(edges <= R_fine(end));
nRegions = length(valid_edges) - 1;
region_colors = {[0.8 0 0.8], [1 0 0], [0 0.7 0], [0 0 1], [0.5 0.5 0.5]};
for k = 1:nRegions
    mid_x = (valid_edges(k) + valid_edges(k+1)) / 2;
    if k <= length(label_order)
        text(mid_x, top_y, label_order{k}, 'FontSize', 11, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Color', region_colors{k});
    end
end
hold off;
set(gca, 'YScale', 'log');
grid on;
xlabel('Horizontal Distance R (m)', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs Distance with AM (Urban, h_{UAV}=%dm, A2G Fading, Monte Carlo)', h_uav), ...
    'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
xlim([0 R_fine(end)]);
ylim(ylims);
set(gca, 'FontSize', 11);
fprintf('\nDone.\n');

%% ==================== Functions ====================
function theta = calc_theta(R, h)
    if R == 0
        theta = 90;
    else
        theta = atan(h/R) * 180/pi;
    end
end
function d = calc_d(R, h)
    d = sqrt(R^2 + h^2);
end
function P_los = calc_Plos(R, h, a, b)
    theta = calc_theta(R, h);
    P_los = 1 / (1 + a*exp(-b*(theta - a)));
end
function PL_fs = calc_PLfs(R, h, fc, c)
    d = calc_d(R, h);
    PL_fs = 20*log10(4*pi*fc*d/c);
end
function PL_los = calc_PL_los(R, h, fc, c, L_los)
    PL_fs = calc_PLfs(R, h, fc, c);
    PL_los = PL_fs + L_los;
end
function PL_nlos = calc_PL_nlos(R, h, fc, c, L_nlos)
    PL_fs = calc_PLfs(R, h, fc, c);
    PL_nlos = PL_fs + L_nlos;
end

