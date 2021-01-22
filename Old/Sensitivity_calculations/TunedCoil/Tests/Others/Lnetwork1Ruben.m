% This test program uses a inverted L-matching network to tune the coil and match it to the frequency of interest
%   |--Rs----Cs_tun1----------|\
%   |                 |       | \
%   Ls              Cs_tun2   |Ampl>
%   |                 |       | /
%   |-------------------------|/
clear all, clc, close all
% Define the characteristics of the coil
Rs = 71.7;
Ls = 57.4e-3;
% Define the frequency of interest
f0 = 1e4;
% Define the noise characteristics of the amplifier
en = 1.3e-9; % V/sqrt(Hz)
in = 0.8e-12; % A/sqrt(Hz)
Rap = 60e6; % Input resistance
% Rap = 1625; % Input resistance
Cap = 2e-12; % Input capacitance
% Calculations acording to 'lect7_match.pdf' of Berkeley
w0 = 2*pi*f0;
RnoiseMatch = en/in; % optimal impedance for noise matching
% RnoiseMatch = 160;
Xs = 1j*w0*Ls;
%% For noise matching
m = RnoiseMatch/Rs; % Calculate the boosting factor
% Rpl = RnoiseMatch*Rap/(RnoiseMatch+Rap); % R- Trying to tune for low input impedance amplifiers
% m = Rpl/Rs; % R- Calculate the boosting factor
Qreq = sqrt(m-1); % Compute the required circuit Q by (1 + Q2) = m => Q = sqrt(m-1)
Xs_total = 1j*Qreq*Rs; % Qreq = Xs_t/Rs;
Xs_tun1 = Xs_total - Xs;
Xs_total_p = Xs_total*(1+Qreq^-2); % Equivalent parallel impedance
Xs_tun2 = -Xs_total_p;

% %% For noise fine matching 2
% m = RnoiseMatch/Rs; % Calculate the boosting factor
% Rpl = Rs*Rap/(Rs+Rap); % -R
% Qreq = sqrt(m-1); % Compute the required circuit Q by (1 + Q2) = m => Q = sqrt(m-1)
% Xs_tun2 = -1j*(Rpl/Qreq);
% Qt2_Rap = Rap/abs(Xs_tun2);
% Xs_tun2_s = Xs_tun2*((1+Qt2_Rap^-2));
% Xs_tun1 = -1*(Xs_tun2_s + Xs);
% Rs_p = (Rs*abs(Xs_tun2_s^2))/(Rs^2+abs(Xs_tun2_s^2));

% %% For noise matching 3
% Qreq = sqrt(RnoiseMatch/Rs-1);
% Xs_tun1_temp = 1j*Qreq*Rs; % inductive
% Xs_tun1 = Xs_tun1_temp - Xs;
% Xs_tun2 = -1j*(RnoiseMatch/Qreq);

% Calculate the value of the capacitors needed
Ct1 = 1/(1j*w0*Xs_tun1);
Ct2_temp = 1/(1j*w0*Xs_tun2);
Ct2 = Ct2_temp-Cap; % Consider the input capacitance of the amplifier

disp(['The capacitiors calculated to noise match the amplifier (' num2str(RnoiseMatch) 'ohm) are:']) 
disp(['Ct1 = ' num2str(Ct1) 'f and Ct2 = ' num2str(Ct2) 'f'])

% Visualize the impedance seen from the amplifier (Zth)
fname='L_inv_Zth_from_amp.cir';
scam
Zth_from_amp = -V1/I_V1;
Vref =1; % A reference voltage for the transfer function
[nsym,dsym]=numden(subs(Zth_from_amp, [V1,L1,R1,C1,C2], [Vref,Ls,Rs,Ct1,Ct2]));
% convert them to MATLAB polynomials
mySys=tf(sym2poly(nsym),sym2poly(dsym));
% Bode Plot
figure(2)
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.MagUnits = 'abs';
bodefreq = logspace(log10(2*pi*1e3),log10(2*pi*1e5),10000);
h = bodeplot(mySys,bodefreq,opts);
Zth_w0 = evalfr(mySys,0+w0*1j);
temp1 = 0+w0*1j;
grid on
disp(['The noise match of the amplifier (' num2str(RnoiseMatch) 'ohm) has been matched to ' num2str(Zth_w0)]) 
% Visualize the voltage seen by the amplifier
fname='L_inv.cir';
scam
Vin_amp = v_4;
title('Zth seen from the amplifier')
[nsym2,dsym2]=numden(subs(Vin_amp, [V1,Ls_s,Rs_s,Ct1_s,Ct2_s,Ra_s,Ca_s], [Vref,Ls,Rs,Ct1,Ct2,Rap,Cap]));
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

%% Calculate SNR for different matching impedances (real)
% Visualize for different matching impedances
Rmatch_sweep = logspace(log10(Rs),5,40);
for k = 1:length(Rmatch_sweep)
%   amplifier
    m_sweep(k) = Rmatch_sweep(k)/Rs; % Calculate the boosting factor
    Q_sweep(k) = sqrt(m_sweep(k)-1); % Compute the required circuit Q by (1 + Q2) = m => Q = sqrt(m-1)
    Xs_total_sweep = 1j*Q_sweep(k)*Rs; % Qreq = Xs_t/Rs;
    Xs_tun1_sweep = Xs_total_sweep - Xs; 
    Xs_total_p_sweep = Xs_total_sweep*(1+Q_sweep(k)^-2); % Equivalent parallel impedance
    Xs_tun2_temp_sweep = -Xs_total_p_sweep;
%     Calculate the value of the capacitors needed
    Ct1_sweep(k) = 1/(1j*w0*Xs_tun1_sweep);
    Ct2_temp_sweep = 1/(1j*w0*Xs_tun2_temp_sweep);
    Ct2_sweep(k) = Ct2_temp_sweep-Cap; % Consider the input capacitance of the amplifier
%   Calculate thevening impedance
    [nsym_zth_sweep,dsym_zth_sweep]=numden(subs(Zth_from_amp, [V1,L1,R1,C1,C2], [Vref,Ls,Rs,Ct1_sweep(k),Ct2_sweep(k)]));
    % convert them to MATLAB polynomials
    mySys_zth_sweep=tf(sym2poly(nsym_zth_sweep),sym2poly(dsym_zth_sweep));
    Zth_w0_zth_sweep(k) = evalfr(mySys_zth_sweep,0+w0*1j);
%   Calculate the voltage gain due to the maching network seen at the
%   amplifier
    [nsym_sweep,dsym_sweep]=numden(subs(Vin_amp, [V1,Ls_s,Rs_s,Ct1_s,Ct2_s,Ra_s,Ca_s], [Vref,Ls,Rs,Ct1_sweep(k),Ct2_sweep(k),Rap,Cap]));
%     convert them to MATLAB polynomials
    mySys_sweep=tf(sym2poly(nsym_sweep),sym2poly(dsym_sweep));
    Vin_amp_w0_sweep(k) = evalfr(mySys_sweep,0+w0*1j);    
end
% Calculate the input inpedance of the amplifier
Rtotal = 1./(1./real(Zth_w0_zth_sweep)+1/Rap);
Qplot = abs(Rtotal ./ Xs_tun2_temp_sweep); 
SNR =sqrt(1+Q_sweep.^2)./sqrt(en^2+(in*Rs.*m_sweep).^2); % Maybe we can find the maximum on this simple equation?
figure(4)
semilogx(Q_sweep, SNR);
grid on
hold on
plot(Qreq,abs((Vin_amp_w0/Vref)./sqrt(en^2+(in*Zth_w0).^2)),'*r')
figure(4)
plot(Q_sweep,abs((Vin_amp_w0_sweep./Vref)./sqrt(en^2+(in.*Zth_w0_zth_sweep).^2)),'.b')
xlabel('Q'), ylabel('SNR')
figure(5)
semilogx(Rmatch_sweep,abs((Vin_amp_w0_sweep/Vref)./sqrt(en^2+(in.*Zth_w0_zth_sweep).^2)),'.b')
xlabel('R match'), ylabel('SNR'), grid on, hold on
plot(RnoiseMatch,abs((Vin_amp_w0./Vref)./sqrt(en^2+(in.*Zth_w0).^2)),'*r')