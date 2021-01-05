% This file simulates the awgn noise 

clc; 
clear;
close all; 

% Data construction 
ts = 0.001; 
tb = 1;
snrdb = 10;
t = 0 : ts : tb;
s = cos(2 * pi * t);
ps = sum(s.^2) / length(s);

% Add noise with awgn 
r1 = awgn(s, snrdb, 10 * log10(ps));
n1 = s - r1;
pn1 = sum(n1.^2) / length(n1);

% Add noise explicitly. 
snr = 10^(snrdb / 10);
pn2 = ps / snr;
n2 = sqrt(pn2) * randn(1, length(s));
r2 = s + n2;

% Report 
fprintf('snr1 = %f dB\n', 10 * log10(ps / pn1));
fprintf('snr2 = %f dB\n', 10 * log10(ps / pn2));
fprintf('var1 = %f\n', var(n1));
fprintf('var2 = %f\n', var(n2));

% Plots 
plot(t, s)
hold on 
plot(t, r1)  
hold on 
plot(t, r2)

