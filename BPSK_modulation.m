clc;
clear;
close all;

modOrder = 2;
num_bits = 200000;
k = log2(modOrder);
num_symbols = num_bits / k;
bit_rate = 1e3;
bit_duration = 1 / bit_rate;
symbol_duration = bit_duration * k;
fc = 1e3;
fs = 1e5;
ts = 1/fs;
samples_per_bit = bit_duration/ts;
samples_per_symbol = samples_per_bit * k;

source_bit = randi([0,1],num_bits,1);
source_bit = source_bit.';
fprintf('Binary Data Sequence:\n');
fprintf('%d ', source_bit);

modOut = pskmod(source_bit, modOrder, pi, "gray");
fprintf('\n\nBPSK Symbols:\n');
fprintf('%+d ', modOut);
fprintf('\n\n');

x = (0:modOrder-1);                   % Integer input
symgray = qammod(x,modOrder,'gray','UnitAveragePower',true); % 64-QAM output (Gray-coded)

scatterplot(symgray,1,0,'y*');
for k = 1:modOrder
    text(real(symgray(k)) - 0.0,imag(symgray(k)) + 0.1, dec2base(x(k),2,6),'Color',[0 1 0]);
    text(real(symgray(k)) - 0.1,imag(symgray(k)) + 0.1, num2str(x(k)),'Color',[0 1 0]);
end

xlim([-1.3,1.3]);
ylim([-1.3,1.3]);
title('BPSK Symbol Constellation','FontWeight', 'bold', 'FontSize', 14);
