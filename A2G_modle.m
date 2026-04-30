clc;
clear;
close all;

%% ==================== OFDM Parameters ====================
numCarr = 256;
cycPrefLen = 32;
nullIdx = [1:32, 129, 225:256]';
numDataCarr = numCarr - length(nullIdx);
ofdmSymbols = 5000;

%% ==================== A2G Channel Parameters ====================
fc = 2.4e9;          % Carrier frequency (Hz)
c = 3e8;             % Speed of light (m/s)
h_uav = 100;         % UAV altitude (m)
Ptx_dBm = 20;        % Transmit power (dBm)
B = numCarr * 15e3;  % Bandwidth (Hz)
NF_dB = 5;           % Noise figure (dB)
P_N = -174 + 10*log10(B) + NF_dB;  % Noise power (dBm)

a = 9.61; b = 0.16;  % Environment parameters for P_LOS
L_los = 1;           % LOS additional loss (dB)
L_nlos = 20;         % NLOS additional loss (dB)
K_rician = 10;       % Rician K-factor

R_vec = 100:50:1000; % Horizontal distance vector (m)

%% ==================== A2G Channel Model ====================

% Fig 1: Elevation Angle vs Distance
figure('Position',[50 50 1000 600]);
subplot(2,2,1);
theta_vec = zeros(size(R_vec));
for i = 1:length(R_vec)
    theta_vec(i) = calc_theta(R_vec(i), h_uav);
end
plot(R_vec/1000, theta_vec, 'b-', 'LineWidth', 2);
grid on; xlabel('Horizontal Distance R (km)'); ylabel('Elevation Angle \theta (deg)');
title('Elevation Angle vs Distance');

% Fig 2: P_LOS vs Distance
subplot(2,2,2);
Plos_vec = zeros(size(R_vec));
for i = 1:length(R_vec)
    Plos_vec(i) = calc_Plos(R_vec(i), h_uav, a, b);
end
plot(R_vec/1000, Plos_vec, 'r-', 'LineWidth', 2);
grid on; xlabel('Horizontal Distance R (km)'); ylabel('P_{LOS}');
title('LOS Probability vs Distance');
ylim([0 1]);

% Fig 3: Path Loss vs Distance (LOS/NLOS/Avg)
subplot(2,2,3);
PL_los_vec = zeros(size(R_vec));
PL_nlos_vec = zeros(size(R_vec));
PL_avg_vec = zeros(size(R_vec));
for i = 1:length(R_vec)
    PL_los_vec(i) = calc_PL_los(R_vec(i), h_uav, fc, c, L_los);
    PL_nlos_vec(i) = calc_PL_nlos(R_vec(i), h_uav, fc, c, L_nlos);
    PL_avg_vec(i) = calc_PL_avg(R_vec(i), h_uav, fc, c, a, b, L_los, L_nlos);
end
plot(R_vec/1000, PL_los_vec, 'b-', 'LineWidth', 2); hold on;
plot(R_vec/1000, PL_nlos_vec, 'r-', 'LineWidth', 2);
plot(R_vec/1000, PL_avg_vec, 'k--', 'LineWidth', 2);
grid on; xlabel('Horizontal Distance R (km)'); ylabel('Path Loss (dB)');
legend('PL_{LOS}','PL_{NLOS}','PL_{avg}','Location','southeast');
title('Path Loss vs Distance');

% Fig 4: SNR vs Distance (LOS/NLOS/Avg)
subplot(2,2,4);
SNR_los_vec = zeros(size(R_vec));
SNR_nlos_vec = zeros(size(R_vec));
SNR_avg_vec = zeros(size(R_vec));
for i = 1:length(R_vec)
    SNR_los_vec(i) = Ptx_dBm - PL_los_vec(i) - P_N;
    SNR_nlos_vec(i) = Ptx_dBm - PL_nlos_vec(i) - P_N;
    SNR_avg_vec(i) = Ptx_dBm - PL_avg_vec(i) - P_N;
end
plot(R_vec/1000, SNR_los_vec, 'b-', 'LineWidth', 2); hold on;
plot(R_vec/1000, SNR_nlos_vec, 'r-', 'LineWidth', 2);
plot(R_vec/1000, SNR_avg_vec, 'k--', 'LineWidth', 2);
grid on; xlabel('Horizontal Distance R (km)'); ylabel('SNR (dB)');
legend('SNR_{LOS}','SNR_{NLOS}','SNR_{avg}','Location','northeast');
title('Received SNR vs Distance');

%% ==================== Functions ====================

function theta = calc_theta(R, h)
    theta = atan(h/R) * 180/pi;
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

function PL_avg = calc_PL_avg(R, h, fc, c, a, b, L_los, L_nlos)
    P_los = calc_Plos(R, h, a, b);
    PL_los = calc_PL_los(R, h, fc, c, L_los);
    PL_nlos = calc_PL_nlos(R, h, fc, c, L_nlos);
    PL_avg_linear = P_los * 10^(PL_los/10) + (1-P_los) * 10^(PL_nlos/10);
    PL_avg = 10*log10(PL_avg_linear);
end