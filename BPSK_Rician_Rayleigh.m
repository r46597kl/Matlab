clc;
clear;
close all;

modOrder = 2;
k = log2(modOrder);

numCarr = 256;
cycPrefLen = 32;
numGBCarr = numCarr / 8;
gbLeft = 1:numGBCarr;
gbRight = (numCarr-numGBCarr+1):numCarr;
dcIdx = (numCarr/2)+1;
nullIdx = [gbLeft dcIdx gbRight]';
numDataCarr = numCarr - length(nullIdx);

num_bps_ofdm = numDataCarr * k;
num_ofdm_symbols = 50000;
numBits = num_bps_ofdm * num_ofdm_symbols;

srcBits = randi([0,1],numBits,1);
pskmodOut = pskmod(srcBits,modOrder,"InputType","bit");
pskmodOut = reshape(pskmodOut,numDataCarr,[]);

EsN0_dB = 0:1:20;
EsN0 = 10.^(EsN0_dB/10);

K_dB = 10;
K = 10^(K_dB/10);

% For BPSK: Eb/N0 = Es/N0 / k, and k=1, so Eb/N0 = Es/N0
EbN0_dB = EsN0_dB - 10*log10(k);

BPSK_BER_AWGN_theoretical = berawgn(EbN0_dB,'psk',2,'nondiff');
BPSK_BER_Rayleigh_theoretical = berfading(EbN0_dB,'psk',2,1);
BPSK_BER_Rician_theoretical = berfading(EbN0_dB,'psk',2,1,K);

BPSK_BER_AWGN_simulated = zeros(size(EsN0_dB));
BPSK_BER_Rayleigh_simulated = zeros(size(EsN0_dB));
BPSK_BER_Rician_simulated = zeros(size(EsN0_dB));

fprintf('-----OFDM BPSK BER Simulation-----\n');
fprintf('Rician K Factor: %d dB\n\n', K_dB);

%% AWGN Simulation
fprintf('--- AWGN Channel ---\n');
for i = 1:length(EsN0_dB)

    snr_dB = EsN0_dB(i) + 10*log10(numDataCarr/numCarr);

    ofdmModOut = ofdmmod(pskmodOut,numCarr,cycPrefLen,nullIdx);

    rxSignal = awgn(ofdmModOut,snr_dB,"measured");
    ofdmDemodOut = ofdmdemod(rxSignal,numCarr,cycPrefLen,0,nullIdx);

    ofdmDemodOut = ofdmDemodOut(:);
    demodOut = pskdemod(ofdmDemodOut, modOrder, "OutputType","bit");

    isBitError = srcBits ~= demodOut;
    numBitError = nnz(isBitError);
    BPSK_BER_AWGN_simulated(i) = numBitError/numBits;

    fprintf('SNR = %.2f dB\tEs/N0 = %d dB\tBER = %.5f\tError Number = %d\n',snr_dB,EsN0_dB(i),BPSK_BER_AWGN_simulated(i),numBitError);
end

%% Rician Simulation
fprintf('\n--- Rician Channel (K=%d dB) ---\n', K_dB);
for i = 1:length(EsN0_dB)

    snr_dB = EsN0_dB(i) + 10*log10(numDataCarr/numCarr);

    [numDataCarr, num_ofdm_symbols] = size(pskmodOut);

    LOS = sqrt(K / (K + 1));
    scatter = sqrt(1 / (K + 1)) * (randn(numDataCarr,num_ofdm_symbols)+1j*randn(numDataCarr,num_ofdm_symbols))/sqrt(2);
    h_rician = LOS + scatter;

    faded_symbols = h_rician .* pskmodOut;

    ofdmModOut = ofdmmod(faded_symbols,numCarr,cycPrefLen,nullIdx);

    rxSignal = awgn(ofdmModOut,snr_dB,"measured");
    ofdmDemodOut = ofdmdemod(rxSignal,numCarr,cycPrefLen,0,nullIdx);

    ofdmDemodOut_equalized = ofdmDemodOut ./ h_rician;

    ofdmDemodOut_equalized = ofdmDemodOut_equalized(:);
    demodOut = pskdemod(ofdmDemodOut_equalized, modOrder, "OutputType","bit");

    isBitError = srcBits ~= demodOut;
    numBitError = nnz(isBitError);
    BPSK_BER_Rician_simulated(i) = numBitError/numBits;

    fprintf('SNR = %.2f dB\tEs/N0 = %d dB\tBER = %.5f\tError Number = %d\n',snr_dB,EsN0_dB(i),BPSK_BER_Rician_simulated(i),numBitError);
end

%% Rayleigh Simulation
fprintf('\n--- Rayleigh Channel ---\n');
for i = 1:length(EsN0_dB)

    snr_dB = EsN0_dB(i) + 10*log10(numDataCarr/numCarr);

    [numDataCarr, num_ofdm_symbols] = size(pskmodOut);
    h_rayleigh = (randn(numDataCarr,num_ofdm_symbols)+1j*randn(numDataCarr,num_ofdm_symbols))/sqrt(2);

    faded_symbols = h_rayleigh .* pskmodOut;

    ofdmModOut = ofdmmod(faded_symbols,numCarr,cycPrefLen,nullIdx);

    rxSignal = awgn(ofdmModOut,snr_dB,"measured");
    ofdmDemodOut = ofdmdemod(rxSignal,numCarr,cycPrefLen,0,nullIdx);

    ofdmDemodOut_equalized = ofdmDemodOut ./ h_rayleigh;

    ofdmDemodOut_equalized = ofdmDemodOut_equalized(:);
    demodOut = pskdemod(ofdmDemodOut_equalized, modOrder, "OutputType","bit");

    isBitError = srcBits ~= demodOut;
    numBitError = nnz(isBitError);
    BPSK_BER_Rayleigh_simulated(i) = numBitError/numBits;

    fprintf('SNR = %.2f dB\tEs/N0 = %d dB\tBER = %.5f\tError Number = %d\n',snr_dB,EsN0_dB(i),BPSK_BER_Rayleigh_simulated(i),numBitError);
end

%% Figure 1: Theoretical BER
figure(1);
set(gcf, 'Position', [0 100 950 600]);
semilogy(EsN0_dB,BPSK_BER_AWGN_theoretical,'k--','LineWidth',2,'DisplayName','Theoretical AWGN BER');
hold on;
semilogy(EsN0_dB,BPSK_BER_Rayleigh_theoretical,'r-','LineWidth',2,'DisplayName','Theoretical Rayleigh BER');
semilogy(EsN0_dB,BPSK_BER_Rician_theoretical,'b-','LineWidth',2,'DisplayName',['Theoretical Rician BER (K=' num2str(K_dB) ' dB)']);
hold off;
xlabel('E_s/N_0 (dB)','FontWeight','bold','FontSize',14);
ylabel('Bit Error Rate (BER)','FontWeight','bold','FontSize',14);
title('OFDM BPSK BER Performance Comparisons','FontWeight','bold','FontSize',16);
legend('show','Location','southwest');
xlim([min(EsN0_dB),max(EsN0_dB)]);
ylim([1e-7,1]);
set(gca,'YScale','log');
grid on;

%% Figure 2: Simulated BER
figure(2);
set(gcf, 'Position', [975 100 950 600]);
semilogy(EsN0_dB,BPSK_BER_AWGN_simulated,'ko-','LineWidth',2,'MarkerFaceColor','black','DisplayName','Simulated AWGN BER');
hold on;
semilogy(EsN0_dB,BPSK_BER_Rayleigh_simulated,'ro-','LineWidth',2,'MarkerFaceColor','red','DisplayName','Simulated Rayleigh BER');
semilogy(EsN0_dB,BPSK_BER_Rician_simulated,'bo-','LineWidth',2,'MarkerFaceColor','blue','DisplayName',['Simulated Rician BER (K=' num2str(K_dB) ' dB)']);
hold off;
xlabel('E_s/N_0 (dB)','FontWeight','bold','FontSize',14);
ylabel('Bit Error Rate (BER)','FontWeight','bold','FontSize',14);
title('OFDM BPSK BER Performance Comparisons','FontWeight','bold','FontSize',16);
legend('show','Location','southwest');
xlim([min(EsN0_dB),max(EsN0_dB)]);
ylim([1e-7,1]);
set(gca,'YScale','log');
grid on;