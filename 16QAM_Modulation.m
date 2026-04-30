clc;
clear;
close all;

modOrder = 16;
num_bits = 10000;
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

display_bits = 20;
display_symbols = display_bits/k;
display_time = display_symbols*symbol_duration;

source_bit = randi([0,1],num_bits,1);
fprintf('Binary Data Sequence:\n');
fprintf('%d ', source_bit);

modOut = qammod(source_bit,modOrder,'InputType','bit','UnitAveragePower',true);

bit_groups = reshape(source_bit,k, num_symbols)';
I_bits = bit_groups(:,1:2);
Q_bits = bit_groups(:,3:4);

symbol = bi2de(bit_groups, 'left-msb');

total_time = num_symbols*symbol_duration;
t_total = ts:ts:total_time;

carrier_I = cos(2*pi*fc*t_total);
carrier_Q = sin(2*pi*fc*t_total);

carrier_display_time = 10;
carrier_samples = carrier_display_time*fs/1000;

qam16_I = [];
qam16_Q = [];
   
for i = 1:num_symbols
    I = real(modOut(i));
    Q = imag(modOut(i));
    qam16_I = [qam16_I, I * ones(1,samples_per_symbol)];
    qam16_Q = [qam16_Q, Q * ones(1,samples_per_symbol)];
end

I_dec = bi2de(I_bits,'left-msb');
Q_dec = bi2de(Q_bits,'left-msb');

QAM16_I = qam16_I .* carrier_I;
QAM16_Q = qam16_Q .* carrier_Q;
QAM16_sig = QAM16_I + QAM16_Q;

x = (0:15);                   % Integer input
symgray = qammod(x,modOrder,'gray'); % 16-QAM output (Gray-coded)

scatterplot(symgray,1,0,'y*');
for k = 1:modOrder
    text(real(symgray(k)) - 0.0,imag(symgray(k)) + 0.3, dec2base(x(k),2,4),'Color',[0 1 0]);
    text(real(symgray(k)) - 0.5,imag(symgray(k)) + 0.3, num2str(x(k)),'Color',[0 1 0]);
end

xlim([-4,4]);
ylim([-4,4]);
title('16-QAM Symbol Constellation','FontWeight','bold','FontSize',14);
