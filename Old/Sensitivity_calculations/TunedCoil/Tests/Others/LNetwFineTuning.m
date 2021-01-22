%% Accurate capacitor calculation
clear all, close all, clc
Rs = 71.7;
Ls = 57.4e-3;
Ra = 60e6;
Ca = 2e-12;
% V1_num = 1;
f0 = 1e4;
w0 = 2*pi*f0;
en_amp = 1.3e-9;
in_amp = 0.8e-12;
fname='L_inv_Zth_from_amp_simp.cir'; % This model takes the inductor and the series capacitor as a total inductor
scam
% First condition is that Zth should be = Rmatch
Zth_s = -V1/I_V1;
Rmatch = en_amp/in_amp;
% Second condition is that it should be in resonance
syms wt R1t L1t C2t R2t real
% Xs = 1j*w*L1;
Xs_mag = wt*L1t;
Xs_p = 1j*((R1t^2+Xs_mag^2)/Xs_mag);
Xc = 1/(1j*wt*C2t);
% Redefine 's' in terms of 'jw'
% Zth = subs(Zth_s, [s,R1,L1,C2],[(1j*wt),R1t,L1t,C2t]); % s = 1j*w;
Zth_t = subs(Zth_s, [s,R1,L1,C2],[(1j*w0),Rs,L1t,C2t]); % s = 1j*w;
Xs_p_t = subs(Xs_p,[wt,R1t],[w0,Rs]);
Xc_t = subs(Xc,[wt],[w0]);
% % Convert equations into 'symfun'format
% Zth_f = symfun(Zth,[w,R1,L1,C2]);
% Xs_p_f = symfun(Xs_p,[w,R1,L1]);
% Xc_f = symfun(Xc,[w,C2]);
% Preapare the functions for the solver
% f1 = Xs_p+Xc;
% f2 = Zth-Rmatch;
f1_t = Xs_p_t+Xc_t;
f2_t = Zth_t-Rmatch;
% Find capacitors for the resonant point and Rmatch
% val1 = solve(f1==0, f2==0,'L1t','C2t');
val1_t = solve(f1_t==0, f2_t==0,'L1t','C2t');
L1tot = double(val1_t.L1t(2));
C2tot = double(val1_t.C2t(2));

Xs_p_f = symfun(Xs_p_t,[L1t]);
Xc_f = symfun(Xc_t,[C2t]);
Zth_f = symfun(Zth_t,[L1t,C2t]);

double(Xs_p_f(L1tot))
double(Xc_f(C2tot))
double(Zth_f(L1tot,C2tot))

Rs_p = (71.7^2+333.7239^2)/71.7;

% plot the impedance bode
figure(2);
h = plotL_inv(Zth_s,[R1,L1,C2],[Rs,L1tot,C2tot]);
title('Impedance seen from the amplifier')
% plote the voltage amplification
fname='L_inv_simp.cir';
scam
Vin_amp = v_3;
figure(3)
h = plotL_inv(Vin_amp,[V1,L1,R1,C2,R2], [1,L1tot,Rs,C2tot,Ra]);
title('Voltage seen at the amplifier')
hold on
% Equations to fine tune from Vin_amp
Vin_amp_f = symfun(Vin_amp,[s,V1,R1,L1,C2,R2]);
dVin_amp_f = diff(Vin_amp_f,s);
w0_t3 = solve(dVin_amp_f==0,'s');
dVin_amp = diff(Vin_amp,'s');
w0_t = solve(dVin_amp==0,'s');
f0_comp  = double(subs(w0_t,[R1,L1,C2,R2],[Rs,L1tot,C2tot,Ra]))/(2*pi); % check if it is right
Vin_amp2 = subs(Vin_amp,[s,V1,R1,L1,C2,R2],[1j*wt,1,R1t,L1t,C2t,R2t]);
dVin_amp2 = diff(Vin_amp2,wt);
w0_t2 = solve(dVin_amp2==0,'wt');
f0_comp2  = double(subs(w0_t2,[R1t,L1t,C2t,R2t],[Rs,L1tot,C2tot,Ra]))/(2*pi); % check if it is right

% Fine tuning by eye
% TunL = +1e-4; % Up and back
% TunC = 0;
% TunL = 0; 
% TunC = 0.7e-9; % Down and back
TunL = 2.93e-5; % Up and back
TunC = 7.8e-10; % Down and back
L1tot_t = L1tot+TunL;
C2tot_t = C2tot+TunC;
figure(2);
hold on
h = plotL_inv(Zth_s,[R1,L1,C2],[Rs,L1tot_t,C2tot_t]);
title('Impedance seen from the amplifier')
% plote the voltage amplification
figure(3);
hold on
h = plotL_inv(Vin_amp,[V1,R1,L1,C2,R2], [1,Rs,L1tot_t,C2tot_t,Ra]);
title('Voltage seen at the amplifier')
Xc1 = w0*(Ls-L1tot_t);
C1_val_fine = -1/(w0*Xc1); % C1_val = 1/(w0^2*(L1tot_t-Ls)); % 1/wC = -wL
C2_val_fine = C2tot_t-Ca;
