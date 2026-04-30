clc;
clear;
close all;

%% ==================== Parameters ====================
numCarr = 256;
cycPrefLen = 32;
nullIdx = [1:32, 129, 225:256]';
numDataCarr = numCarr - length(nullIdx);
ofdmSymbols = 5000;

fc = 2.4e9;
c = 3e8;
Ptx_dBm = 30;
B = numCarr * 15e3;
NF_dB = 5;
P_N = -174 + 10*log10(B) + NF_dB;

a = 9.61; b = 0.16;
L_los = 1;
L_nlos = 20;
K_rician = 10;
BER_target = 1e-2;

% BPSK
modOrder = 2;
k_mod = log2(modOrder);
numBitsPerMC = numDataCarr * k_mod * ofdmSymbols;

h_vec = [50, 100, 200, 300];
R_vec = 0:50:1000;
numMC = 30;

colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.929 0.694 0.125], [0.494 0.184 0.556]};
markers = {'o', 's', 'd', '^'};

%% ==================== Simulation ====================
figure('Position',[100 100 900 600]);

for hi = 1:length(h_vec)
    h_uav = h_vec(hi);
    BER_h = zeros(size(R_vec));
    fprintf('--- h_uav = %d m ---\n', h_uav);
    
    for ri = 1:length(R_vec)
        R = R_vec(ri);
        P_los = calc_Plos(R, h_uav, a, b);
        PL_los_val = calc_PL_los(R, h_uav, fc, c, L_los);
        PL_nlos_val = calc_PL_nlos(R, h_uav, fc, c, L_nlos);
        
        totalErrors = 0;
        totalBits = 0;
        
        for mc = 1:numMC
            srcBits = randi([0,1], numBitsPerMC, 1);
            modOut = reshape(pskmod(srcBits, modOrder, 'InputType', 'bit'), numDataCarr, []);
            
            % --- Per-OFDM-symbol LOS/NLOS decision ---
            isLOS = rand(1, ofdmSymbols) < P_los;
            
            h_ch = zeros(numDataCarr, ofdmSymbols);
            PL_per_sym = zeros(1, ofdmSymbols);
            
            % LOS symbols: Rician fading
            nLOS = sum(isLOS);
            if nLOS > 0
                h_ch(:, isLOS) = sqrt(K_rician/(K_rician+1)) + ...
                    sqrt(1/(K_rician+1)/2) * (randn(numDataCarr, nLOS) + 1j*randn(numDataCarr, nLOS));
                PL_per_sym(isLOS) = PL_los_val;
            end
            
            % NLOS symbols: Rayleigh fading
            nNLOS = sum(~isLOS);
            if nNLOS > 0
                h_ch(:, ~isLOS) = (randn(numDataCarr, nNLOS) + 1j*randn(numDataCarr, nNLOS)) / sqrt(2);
                PL_per_sym(~isLOS) = PL_nlos_val;
            end
            
            % OFDM transmission per symbol
            rxBits_all = zeros(numBitsPerMC, 1);
            bitIdx = 0;
            
            for sym = 1:ofdmSymbols
                dataSym = modOut(:, sym);
                h_sym = h_ch(:, sym);
                SNR_dB = Ptx_dBm - PL_per_sym(sym) - P_N;
                SNR_lin = 10^(SNR_dB/10);

                txFaded = h_sym .* dataSym;
                
                % Add noise
                noisePow = mean(abs(txFaded).^2) / SNR_lin;
                noise = sqrt(noisePow/2) * (randn(numDataCarr,1) + 1j*randn(numDataCarr,1));
                rxSig = txFaded + noise;
                
                % equalization
                rxEq = rxSig ./ h_sym;
                
                rxBits_sym = pskdemod(rxEq, modOrder, 'OutputType', 'bit');
                rxBits_all(bitIdx+1 : bitIdx+numDataCarr*k_mod) = rxBits_sym;
                bitIdx = bitIdx + numDataCarr * k_mod;
            end
            
            totalErrors = totalErrors + nnz(srcBits ~= rxBits_all);
            totalBits = totalBits + numBitsPerMC;
        end
        
        BER_h(ri) = max(totalErrors / totalBits, 1e-7);
        fprintf('  R=%4dm  P_LOS=%.3f  BER=%.2e\n', R, P_los, BER_h(ri));
    end
    
    semilogy(R_vec, BER_h, ['-' markers{hi}], 'Color', colors{hi}, ...
        'LineWidth', 2, 'MarkerSize', 6);
    hold on;
end

yline(BER_target, 'k--', 'BER_{target} = 10^{-2}', 'LineWidth', 1.5);
grid on;
xlabel('Horizontal Distance R (m)', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title('BER vs Distance (Urban, Different UAV Heights, BPSK)', 'FontSize', 14);
legend([arrayfun(@(h) sprintf('h_{UAV} = %d m', h), h_vec, 'UniformOutput', false), ...
    {'BER_{target}'}], 'Location', 'southeast', 'FontSize', 11);
xlim([0 1000]);
ylim([1e-7 1e-1]);
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
