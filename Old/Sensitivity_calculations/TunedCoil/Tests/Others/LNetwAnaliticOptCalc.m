%% Analytically calculate maxima in the SNR for inverted L matching network
clear all, close all, clc
%% Observed equation, ONLY applicable for high input impedance amplifiers
% m = 1+Q^2;
% Some values
e_nnum = 1.3e-9;
i_nnum = 0.8e-12;
Rs_num = 71.7;
Q_sweep = logspace(log10(1), log10(100),100);
Rp_sweep_num = (1+Q_sweep.^2).*Rs_num;
syms Q e_n i_n Rs Rm real
SNR =sqrt(1+Q^2)/sqrt(e_n^2+(i_n*Rs*(1+Q^2)).^2);
Qt = sqrt((Rm/Rs)-1);
SNRt = subs(SNR,Q,Qt);
dSNR = diff(SNR,Q);
dSNRt = diff(SNRt,Rm);
M = solve(dSNR==0,Q);
Mt = solve(dSNRt==0,Rm); % This points at that the optimum is to noise match!!!
for k = 1:length(Q_sweep)
    SNRtemp(k) = eval(subs(SNR,[e_n i_n Rs Q],[e_nnum i_nnum Rs_num Q_sweep(k)]));
    dSNRtemp(k) = eval(subs(dSNR,[e_n i_n Rs Q],[e_nnum i_nnum Rs_num Q_sweep(k)]));
end
figure(1)
subplot(2,1,1)
semilogx(Q_sweep,SNRtemp)
ylabel('SNRtemp'), xlabel('Q')
grid on
subplot(2,1,2)
semilogx(Rp_sweep_num,SNRtemp)
grid on
figure(2)
subplot(2,1,1)
semilogx(Q_sweep,dSNRtemp)
ylabel('dSNRtemp'),xlabel('Q')
grid on
subplot(2,1,2)
semilogx(Rp_sweep_num,SNRtemp)
grid on
Qopt = M(2);
Rp = (1+Qopt^2)*Rs;
Qopt_num= eval(subs(Qopt,[e_n i_n Rs],[e_nnum i_nnum Rs_num]));
Rp_num= eval(subs(Rp,[e_n i_n Rs],[e_nnum i_nnum Rs_num]));

%% Calculated equation
syms Rs Ls Ra Ca V1 w0 real
s = 0+w0*1j;
Ra_num = 60e6;
Ca_num = 2e-12;
V1_num = 1;
w0_num = 1e4;
Ls_num = 57.4e-3;
% Xs = 1j*w0*Ls;
Xs = s*Ls;
Xs_total = 1j*Q*Rs; % Qreq = Xs_t/Rs;
Xs_tun1 = Xs_total - Xs; 
Xs_total_p = Xs_total*(1+Q^-2); % Equivalent parallel impedance
Xs_tun2_temp = -Xs_total_p;
% Calculate the value of the capacitors needed
% Ct1 = 1/(1j*w0*Xs_tun1);
Ct1 = 1/(s*Xs_tun1);
% Ct2_temp = 1/(1j*w0*Xs_tun2_temp);
Ct2_temp = 1/(s*Xs_tun2_temp);
Ct2 = Ct2_temp-Ca; % Consider the input capacitance of the amplifier
Vin_amp = -(Ct1*Ra*V1*s)/(Ca*Ra*s + Ct1*Ra*s + Ct2*Ra*s + Ct1*Rs*s + Ct1*Ls*s^2 + Ca*Ct1*Ls*Ra*s^3 + Ct1*Ct2*Ls*Ra*s^3 + Ca*Ct1*Ra*Rs*s^2 + Ct1*Ct2*Ra*Rs*s^2 + 1);
Zth = (Ct1*Ls*s^2 + Ct1*Rs*s + 1)/(s*(Ct1*Ct2*Ls*s^2 + Ct1*Ct2*Rs*s + Ct1 + Ct2));
SNR2 = simplify((Vin_amp/V1)./sqrt(e_n^2+(i_n*Zth).^2));
dSNR2 = diff(SNR2,Q);
M = solve(dSNR2==0,Q);
for k = 1:length(Q_sweep)
    SNRtemp2(k) = eval(subs(SNR2,[e_n i_n Rs Ls Ra Ca V1 w0 Q],[e_nnum i_nnum Rs_num Ls_num Ra_num Ca_num V1_num w0_num Q_sweep(k)]));
    dSNRtemp2(k) = eval(subs(dSNR2,[e_n i_n Rs Ls Ra Ca V1 w0 Q],[e_nnum i_nnum Rs_num Ls_num Ra_num Ca_num V1_num w0_num Q_sweep(k)]));
end
figure(3)
subplot(2,1,1)
semilogx(Q_sweep,abs(SNRtemp2))
ylabel('SNRtemp'), xlabel('Q')
grid on
subplot(2,1,2)
semilogx(Rp_sweep_num,abs(SNRtemp2))
grid on

%% Similar to the initial theory of SNR but including the real Q of the system accounting for the input impedance of the pre-amplifier
% ---Rsp-------------
%          |    |    |
%          Lp   Cp   Rl
%          |    |    |
% --------------------
clear all,
syms Q e_n i_n Rs Rm Ramp real
% Ramp = 60; Rs = 71.7; e_n = 1.3e-9; i_n = 0.8e-12; Rm = e_n/i_n;
Qmatch = sqrt((Rm/Rs)-1);
% Rsp = (1+Qmatch^2)*Rs;
% Rsp = (1+Rm/Rs-1)*Rs; => Rsp = Rm;
Rsp = Rm;
Xcp = Rsp/Qmatch; % Equivalent parallel impedance
Xlp = Xcp;
Rptot = Rsp*Ramp/(Rsp+Ramp);
Qtotal = Rptot/Xcp;
% Find optimum 'Qmatch'
SNR =sqrt(1+Qtotal^2)/sqrt(e_n^2+(i_n*Rm)^2);
dSNR = simplify(diff(SNR,Rm));
% Mt = solve(dSNR==0,Rm); % This points at that the optimum is to noise match!!!
% Visualize
Ramp_n = 60;
Rs_n = 71.7; 
e_n_n = 1.3e-9;
i_n_n = 0.8e-12;
Rmatch_sweep = logspace(log10(Rs_n),5,40);
%     [nsym_sweep,dsym_sweep]=numden(subs(SNR, [Ramp, Rm, Rs, i_n, e_n], [Ramp_n, Rmatch_sweep(k), Rs_n, i_n_n, e_n_n]));
SNRt = subs(SNR, [Ramp, Rs, i_n, e_n], [Ramp_n, Rs_n, i_n_n, e_n_n]);
SNRf = simplify(symfun(SNRt,[Rm]));
SNRvisu = double(SNRf(Rmatch_sweep));
figure(12)
semilogx(Rmatch_sweep,SNRvisu)
xlabel('R match'), grid on


