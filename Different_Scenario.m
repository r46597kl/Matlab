clc;
clear;
close all;

%% Parameters
numCarr = 256;
cycPrefLen = 32;
nullIdx = [1:32, 129, 225:256]';
numDataCarr = numCarr - length(nullIdx);
ofdmSymbols = 20000;

fc = 2.4e9;
c = 3e8;
h_uav = 300;
Ptx_dBm = 20;
B = numCarr * 15e3;
NF_dB = 5;
P_N = -174 + 10*log10(B) + NF_dB;

K_rician = 10;
BER_target = 1e-2;

%% Environments
envNames  = {'Suburban', 'Urban', 'Dense Urban', 'High-rise Urban'};
a_env     = [4.88,  9.61,  12.08,  27.23];
b_env     = [0.43,  0.16,  0.11,   0.08];
Llos_env  = [0.1,   1.0,   1.6,    2.3];
Lnlos_env = [21,    20,    23,     34];
nScen = length(envNames);

modSchemes = struct( ...
    'name',  {'BPSK', 'QPSK', '16-QAM', '64-QAM'}, ...
    'order', {2,      4,      16,       64}, ...
    'k',     {1,      2,      4,        6}, ...
    'type',  {'psk',  'psk',  'qam',    'qam'} ...
);
nMod = length(modSchemes);

EsN0_dB_vec = 0:2:40;
numSymPerTrial = 2000;
maxTrials = 100;

R_ref = 500;
a_ref = 9.61; b_ref = 0.16;
P_los_ref = calc_Plos(R_ref, h_uav, a_ref, b_ref);

ber_esn0 = zeros(nMod, length(EsN0_dB_vec));

for s = 1:nMod
    M       = modSchemes(s).order;
    bps     = modSchemes(s).k;
    modType = modSchemes(s).type;

    for i = 1:length(EsN0_dB_vec)
        EsN0_dB  = EsN0_dB_vec(i);
        EsN0_lin = 10^(EsN0_dB / 10);
        totalErrors = 0;
        totalBits   = 0;
        numBits_trial = numSymPerTrial * bps;

        for trial = 1:maxTrials
            srcBits = randi([0, 1], numBits_trial, 1);

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

            rxEq = rxSig ./ h;

            if strcmp(modType, 'psk')
                rxBits = pskdemod(rxEq, M, 'OutputType', 'bit');
            else
                rxBits = qamdemod(rxEq, M, 'OutputType', 'bit', 'UnitAveragePower', true);
            end

            totalErrors = totalErrors + nnz(srcBits ~= rxBits);
            totalBits   = totalBits + numBits_trial;
        end

        ber_esn0(s, i) = max(totalErrors / totalBits, 1e-8);
        fprintf('Es/N0=%2ddB  Trials=%2d  Bits=%.0e  Errors=%6d  BER=%.2e\n', ...
            EsN0_dB, maxTrials, totalBits, totalErrors, ber_esn0(s, i));
    end
end

am_thresh = zeros(1, nMod);
for s = 1:nMod
    ber_s = ber_esn0(s, :);
    valid = ber_s > 0 & ber_s < 0.5;
    [x_u, ia] = unique(log10(ber_s(valid)), 'stable');
    y_u = EsN0_dB_vec(valid);
    y_u = y_u(ia);
    am_thresh(s) = interp1(x_u, y_u, log10(BER_target), 'pchip');
end

for s = 1:nMod
    fprintf('  %-7s : %.2f dB\n', modSchemes(s).name, am_thresh(s));
end

R_mc  = 0:50:1000;
nDist = length(R_mc);
numMC_los  = 50;
numMC_nlos = 50;

colors  = {[0.466 0.674 0.188], [0 0.447 0.741], ...
           [0.85 0.325 0.098],  [0.494 0.184 0.556]};
markers = {'s', 'o', 'd', '^'};

BER_mod = zeros(nMod, nScen, nDist);

for ei = 1:nScen
    ai = a_env(ei);
    bi = b_env(ei);
    L_los_i  = Llos_env(ei);
    L_nlos_i = Lnlos_env(ei);

    tic;

    for ri = 1:nDist
        R = R_mc(ri);
        P_los       = calc_Plos(R, h_uav, ai, bi);
        PL_los_val  = calc_PL_los(R, h_uav, fc, c, L_los_i);
        PL_nlos_val = calc_PL_nlos(R, h_uav, fc, c, L_nlos_i);

        SNR_dB_los  = Ptx_dBm - PL_los_val  - P_N;
        SNR_dB_nlos = Ptx_dBm - PL_nlos_val - P_N;

        for s = 1:nMod
            M       = modSchemes(s).order;
            bps     = modSchemes(s).k;
            modType = modSchemes(s).type;
            nBits   = numDataCarr * bps * ofdmSymbols;

            % LOS
            errors_los = 0;
            bits_los   = 0;
            for mc = 1:numMC_los
                h_ch = sqrt(K_rician/(K_rician+1)) + ...
                    sqrt(1/(K_rician+1)/2) * ...
                    (randn(numDataCarr, ofdmSymbols) + 1j*randn(numDataCarr, ofdmSymbols));

                srcBits = randi([0, 1], nBits, 1);

                if strcmp(modType, 'psk')
                    modSym = pskmod(srcBits, M, 'InputType', 'bit');
                else
                    modSym = qammod(srcBits, M, 'InputType', 'bit', ...
                             'UnitAveragePower', true);
                end
                modSym = reshape(modSym, numDataCarr, ofdmSymbols);

                txFaded = ofdmmod(h_ch .* modSym, numCarr, cycPrefLen, nullIdx);
                rxSig   = awgn(txFaded, SNR_dB_los, 'measured');
                rxSym   = ofdmdemod(rxSig, numCarr, cycPrefLen, 0, nullIdx) ./ h_ch;

                if strcmp(modType, 'psk')
                    rxBits = pskdemod(rxSym, M, 'OutputType', 'bit');
                else
                    rxBits = qamdemod(rxSym, M, 'OutputType', 'bit', ...
                             'UnitAveragePower', true);
                end

                errors_los = errors_los + nnz(srcBits ~= rxBits(:));
                bits_los   = bits_los + nBits;
            end
            ber_los = errors_los / bits_los;

            % NLOS
            errors_nlos = 0;
            bits_nlos   = 0;
            for mc = 1:numMC_nlos
                h_ch = (randn(numDataCarr, ofdmSymbols) + ...
                    1j*randn(numDataCarr, ofdmSymbols)) / sqrt(2);

                srcBits = randi([0, 1], nBits, 1);

                if strcmp(modType, 'psk')
                    modSym = pskmod(srcBits, M, 'InputType', 'bit');
                else
                    modSym = qammod(srcBits, M, 'InputType', 'bit', ...
                             'UnitAveragePower', true);
                end
                modSym = reshape(modSym, numDataCarr, ofdmSymbols);

                txFaded = ofdmmod(h_ch .* modSym, numCarr, cycPrefLen, nullIdx);
                rxSig   = awgn(txFaded, SNR_dB_nlos, 'measured');
                rxSym   = ofdmdemod(rxSig, numCarr, cycPrefLen, 0, nullIdx) ./ h_ch;

                if strcmp(modType, 'psk')
                    rxBits = pskdemod(rxSym, M, 'OutputType', 'bit');
                else
                    rxBits = qamdemod(rxSym, M, 'OutputType', 'bit', ...
                             'UnitAveragePower', true);
                end

                errors_nlos = errors_nlos + nnz(srcBits ~= rxBits(:));
                bits_nlos   = bits_nlos + nBits;
            end
            ber_nlos = errors_nlos / bits_nlos;

            BER_mod(s, ei, ri) = max(P_los * ber_los + (1 - P_los) * ber_nlos, 1e-7);
        end

        fprintf('  R=%4dm  P_LOS=%.3f  | BPSK=%.2e  QPSK=%.2e  16QAM=%.2e  64QAM=%.2e\n', ...
            R, P_los, BER_mod(1,ei,ri), BER_mod(2,ei,ri), ...
            BER_mod(3,ei,ri), BER_mod(4,ei,ri));
    end

end

%% AM curve

R_fine = 0:2:R_mc(end);
nFine  = length(R_fine);

BER_smooth_all = zeros(nMod, nScen, nFine);
BER_am_all     = NaN(nScen, nFine);
ModSelect_all  = zeros(nScen, nFine);
SE_all         = zeros(nScen, nFine);
cross_dist_all = zeros(nScen, nMod);

for ei = 1:nScen

    for s = 1:nMod
        ber_raw = squeeze(BER_mod(s, ei, :))';

        for i = 2:nDist
            ber_raw(i) = max(ber_raw(i), ber_raw(i-1));
        end

        BER_smooth_all(s, ei, :) = ...
            10.^interp1(R_mc, log10(ber_raw), R_fine, 'pchip');
    end

    for s = 1:nMod
        ber_s = squeeze(BER_smooth_all(s, ei, :))';
        idx = find(ber_s >= BER_target, 1, 'first');
        if isempty(idx)
            cross_dist_all(ei, s) = R_fine(end) + 1;
        else
            cross_dist_all(ei, s) = R_fine(idx);
        end
    end

    fprintf('\n%s crossing distances:  64QAM=%dm  16QAM=%dm  QPSK=%dm  BPSK=%dm\n', ...
        envNames{ei}, cross_dist_all(ei,4), cross_dist_all(ei,3), ...
        cross_dist_all(ei,2), cross_dist_all(ei,1));

    for i = 1:nFine
        R = R_fine(i);
        selected = 0;
        for s = nMod:-1:1
            if R < cross_dist_all(ei, s)
                selected = s;
                break;
            end
        end

        if selected > 0
            BER_am_all(ei, i)    = BER_smooth_all(selected, ei, i);
            ModSelect_all(ei, i) = selected;
            SE_all(ei, i)        = modSchemes(selected).k;
        end
    end
end

switchLabels = {'64-QAM -> 16-QAM', '16-QAM -> QPSK', 'QPSK -> BPSK', 'BPSK -> No Tx'};
for ei = 1:nScen
    fprintf('\n%s:\n', envNames{ei});
    for s = nMod:-1:1
        cd = cross_dist_all(ei, s);
        if cd <= R_fine(end)
            fprintf('  %s  at R = %d m\n', switchLabels{nMod - s + 1}, cd);
        else
            fprintf('  %s  not reached within %d m\n', switchLabels{nMod - s + 1}, R_fine(end));
        end
    end
end

%% Plot
mod_colors = {[0 0.447 0.741], [0.466 0.674 0.188], ...
              [0.85 0.325 0.098], [0.635 0.078 0.184]};
mod_styles = {'-', '--', '-.', ':'};

figure('Position', [50 50 1200 900]);
for ei = 1:nScen
    subplot(2, 2, ei);
    hold on;

    for s = 1:nMod
        ber_s = squeeze(BER_smooth_all(s, ei, :))';
        plot(R_fine, ber_s, mod_styles{s}, ...
            'Color', mod_colors{s}, 'LineWidth', 1.5, ...
            'DisplayName', modSchemes(s).name);
    end

    ber_am = BER_am_all(ei, :);
    valid = ~isnan(ber_am);
    plot(R_fine(valid), ber_am(valid), '-k', 'LineWidth', 2.5, ...
        'DisplayName', 'AM');

    yline(BER_target, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    for s = 1:nMod
        cd = cross_dist_all(ei, s);
        if cd <= R_fine(end)
            xline(cd, ':', 'Color', mod_colors{s}, 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
    end

    hold off;
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('R (m)');
    ylabel('BER');
    title(envNames{ei}, 'FontSize', 12);
    legend('Location', 'southeast', 'FontSize', 8);
    ylim([1e-6, 1]);
    xlim([0, R_fine(end)]);
end
sgtitle('BER vs Distance with AM', ...
    'FontSize', 14);

%% AM BER
figure('Position', [100 100 900 600]);
hold on;
for ei = 1:nScen
    ber_am = BER_am_all(ei, :);
    valid = ~isnan(ber_am);
    plot(R_fine(valid), ber_am(valid), ['-' markers{ei}], ...
        'Color', colors{ei}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerIndices', 1:50:sum(valid), ...
        'DisplayName', envNames{ei});
end
yline(BER_target, 'k--', 'LineWidth', 1.5, 'DisplayName', 'BER_{target} = 10^{-2}');
hold off;
set(gca, 'YScale', 'log');
grid on;
xlabel('Horizontal Distance R (m)', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title('BER vs Distance with AM', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 11);
xlim([0 R_fine(end)]);
ylim([1e-6 1e-1]);
set(gca, 'FontSize', 11);

%% Modulation choice
figure('Position', [120 80 900 600]);
hold on;
for ei = 1:nScen
    plot(R_fine, ModSelect_all(ei, :), ['-' markers{ei}], ...
        'Color', colors{ei}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerIndices', 1:50:nFine, ...
        'DisplayName', envNames{ei});
end
hold off;
grid on;
xlabel('Horizontal Distance R (m)', 'FontSize', 13);
ylabel('Modulation Selection', 'FontSize', 13);
title('Modulation Selection vs Distance', 'FontSize', 14);
yticks(0:4);
yticklabels({'No Tx', 'BPSK', 'QPSK', '16-QAM', '64-QAM'});
legend('Location', 'northeast', 'FontSize', 11);
xlim([0 R_fine(end)]);
ylim([0 4.5]);
set(gca, 'FontSize', 11);

%% Spectral efficiency
figure('Position', [140 60 900 600]);
hold on;
for ei = 1:nScen
    plot(R_fine, SE_all(ei, :), ['-' markers{ei}], ...
        'Color', colors{ei}, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerIndices', 1:50:nFine, ...
        'DisplayName', envNames{ei});
end
hold off;
grid on;
xlabel('Horizontal Distance R (m)', 'FontSize', 13);
ylabel('Spectral Efficiency (bit/s/Hz)', 'FontSize', 13);
title('Spectral Efficiency vs Distance', ...
    'FontSize', 14);
yticks([0 1 2 4 6]);
yticklabels({'0', '1', '2', '4', '6'});
legend('Location', 'northeast', 'FontSize', 11);
xlim([0 R_fine(end)]);
ylim([0 7]);
set(gca, 'FontSize', 11);

%% Functions

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
