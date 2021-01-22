fname='L_inv_Zth_from_amp_ESR.cir';
% * This circuit allows calculating the Zth seen from the load of an inverted L matching circuit
% * It is done by shifting the excitation to the place of the load
scam
formula1 = V1/I_V1;

fname='L_inv_ESR_Cs.cir';
% * This circuit allows calculating the Zth seen from the load of an inverted L matching circuit
% * It is done by shifting the excitation to the place of the load
scam
formula2 = v_5/V1;


freqs = logspace(1,6,1000);
s_n = freqs*2*pi*1j; 
Rs_n = 407.16e-3;
L_n = 71.003e-6;
Cs_n = 1e-15;
C1_n = 1.0226e-07; 
C2_n = 5.4769e-08;
R_c1_n = 0.11117;
R_c2_n = 0.20756;

fname='L_inv_Zs.cir';
% * This circuit allows calculating the Zth seen from the load of an inverted L matching circuit
% * It is done by shifting the excitation to the place of the load
scam
formula3 = V1/I_V1;
Zth_t = subs(formula3, [Rs,L,Cs],[Rs_n,L_n,Cs_n]); % s = 1j*w;
Zth_f = symfun(Zth_t,s);
Zth_vis = double(Zth_f(s_n));
h = figure;
loglog(freqs,abs(Zth_vis))
grid on, hold on
pos = 100;
plot(abs(s_n(pos))/(2*pi),abs(double(Zth_f(s_n(pos)))),'*r')