%% Analyzing efficiency (power losses)(Q = around 100 in the kHz regime) for an inverted L-network
% Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf.
% See also "Insertion loss". pag 16 of lect7_match.pdf (https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwieq-Tao5fNAhXEW5QKHXKtAEgQFggdMAA&url=http%3A%2F%2Frfic.eecs.berkeley.edu%2F142%2Fpdf%2Flect7_match.pdf&usg=AFQjCNF-SM1zdrBem8kFllXPBTCAM1pQkg)
% Insertion loss = 1-efficiency;
close all, clear all, clc
% Parameters of coil, matching network and amplifier
f = 1e4;
% Coil
Rs = 71.7;
L = 57.4e-3;
% Quality of the matching network
Qc1 = 100;
Qc2 = 100;
% Amplifier
en_amp = 1.3e-9;
in_amp = 0.8e-12;
Ramp = 60e6; % Input resistance
Camp = 2e-12; % Input capacitance

% Calculations
Rmatch = en_amp/in_amp;
Q = sqrt(Rmatch/Rs-1);
% Estimate (rougly) the vaues of the CAPACITORS
[C1, C2, Zth] = L_inv_estim_cap(Rmatch, f, Rs, L, Camp);
% Estimate the power loss
[efficiency2, G_effective2, G_ideal2] = L_inv_effic(Q,Qc1,Qc2);
% Estimate the noise at the input of the amplifier
[en_in_total2, Zth2] = L_inv_noise(en_amp, in_amp, f, Rs, L, C1, Qc1, C2, Qc2, G_effective2);
