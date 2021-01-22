% This test program uses a inverted L-matching network to tune the coil and match it to the frequency of interest
% LCC (inductance-capacitor-capacitor) T network design 
%   |--Rs-------Ct1-------Ct3---|\
%   |                 |         | \
%   Ls              Ct2         |Ampl>
%   |                 |         | /
%   |---------------------------|/
clear all, clc, close all
% Define the characteristics of the coil
Rs = 71.7;
Ls = 57.4e-3;
% Define the frequency of interest abd Bandwidth
f0 = 1e4;
bw = 2e3;
% Define the noise characteristics of the amplifier
en = 1.3e-9; % V/sqrt(Hz)
in = 0.8e-12; % A/sqrt(Hz)
Rap = 60e6; % Input resistance
Cap = 2e-12; % Input capacitance
% Calculate noise matching input inpedance
RnoiseMatch = en/in; % optimal impedance for noise matching
% Calculations acording to
% http://electronicdesign.com/communications/back-basics-impedance-matching-part-3
w0 = 2*pi*f0;
Qreq = f0/bw;
Xl = w0*Ls; % Impedance of the intrinsic inductance of the coil
Xlt = Qreq*Rs; % Inductive impedance needed for the LCC network
Xc1 = Xl-Xlt;  % Capacitor needed to match the intrinsic inductive 
                % imperdance to the matching inductive impedance
Xc3 = RnoiseMatch*sqrt((Rs*(Qreq^2+1)/RnoiseMatch)-1);
Xc2 = (Rs*(Qreq^2+1)/Qreq)*(Qreq*RnoiseMatch/(Qreq*RnoiseMatch-Xc3));
% Calculate the capacitors
Ct1 = 1/(w0*Xc1);
Ct2 = 1/(w0*Xc2);
Ct3 = 1/(w0*Xc3);

% Visualize the impedance seen from the amplifier (Zth)
fname='T_Zth_from_amp.cir';
scam
Zth_from_amp = -V1/I_V1;
Vref =1; % A reference voltage for the transfer function
[nsym,dsym]=numden(subs(Zth_from_amp, [V1,L1,R1,C1,C2,C3], [Vref,Ls,Rs,Ct1,Ct2,Ct3]));
% convert them to MATLAB polynomials
mySys=tf(sym2poly(nsym),sym2poly(dsym));
% Bode Plot
figure(2)
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.MagUnits = 'abs';
bodefreq = logspace(log10(2*pi*1e3),log10(2*pi*1e5),10000);
h = bodeplot(mySys,bodefreq,opts);
title('Zth seen from the amplifier')
Zth_w0 = evalfr(mySys,0+w0*1j);
temp1 = 0+w0*1j;
grid on
disp(['The noise match of the amplifier (' num2str(RnoiseMatch) 'ohm) has been matched to ' num2str(Zth_w0)]) 
% Visualize the voltage seen by the amplifier
fname='T_LCCC.cir';
scam
Vin_amp = v_5;
[nsym2,dsym2]=numden(subs(Vin_amp, [V1,Ls_s,Rs_s,Ct1_s,Ct2_s,Ct3_s,Ra_s,Ca_s], [Vref,Ls,Rs,Ct1,Ct2,Ct3,Rap,Cap]));
% convert them to MATLAB polynomials
mySys2=tf(sym2poly(nsym2),sym2poly(dsym2));
% Bode Plot
figure(3)
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.MagUnits = 'abs';
h = bodeplot(mySys2,bodefreq,opts);
title('Voltage seen at the amplifier')
grid on
Vin_amp_w0 = evalfr(mySys2,0+w0*1j);
disp(['The amplification of the matching network at RnoiseMatch =  (' num2str(RnoiseMatch) 'ohm) is ' num2str(abs(Vin_amp_w0/Vref))]) 
