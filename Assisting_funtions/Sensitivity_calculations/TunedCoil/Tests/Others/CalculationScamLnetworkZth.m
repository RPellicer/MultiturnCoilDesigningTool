clear all, close all, clc
syms w Rs Lm Ct Rm real
Xl = 1j*Lm*w;
Xc = 1/(1j*Ct*w);
Zth_sym = (Xc*(Rs+Xl)/(Xc+(Rs+Xl)));
% convert them to MATLAB polynomials
Zth_f = symfun(Zth_sym,[w,Rs,Lm,Ct]);
Rth_f = simplify(real(Zth_f));
Xth_f = simplify(imag(Zth_f));
dZth_f  = simplify(diff(Zth_f,w));
resonance = solve(dZth_f==0);
reson1 = symfun(resonance(1),[w,Rs,Lm,Ct]);
Rm=Rth_f;
match = solve(Rth_f==Rm,reson1==w, Ct, Lm);
% match = solve(Rs==71.7, reson1==-1e4, Rm==1625, Ct, Lm);
%% Analyzed with SCAM
fname='Zth_seen_from_Rl.cir';
% * This circuit allows calculating the Zth seen from the load of an inverted L matching circuit
% * It is done by shifting the excitation to the place of the load
scam
formula3_temp = V1/I_V1;
formula3 = subs(formula3_temp,[R1,L1,C1],[Rs,Lm,Ct]);
Zth3 = symfun(formula3,[s,Rs,Lm,Ct]);

% % convert the expression to a MATLAB transfer function object
% % First we separate out the numerator and denominator
% % [n,d]=numden(eval(resp));
% % resp = subs(resp,s,w*1i);
% [nsym,dsym]=numden(resp);
% [n,d]=numden(eval(resp));
% % convert them to MATLAB polynomials
% mySys=tf(sym2poly(n),sym2poly(d));

clear all, close all, clc
Rp = 1e11;
Cp = 1e-3;
f = 1e6;
w= 2*pi*f;
Xp = 1/Cp*w;
Rs = (Rp*Xp^2)/((Rp^2)+(Xp^2))
Q = Rp/Xp;
Rs2 = Rp*(1/(1+Q^2))


% % Noise model capacitors
% ESL = 2*pi*f*L;
% ESR = D.F./(2*PI*f*C); %D.F. = dissipation factor
% Z = sqrt((ESR)^2+(ESL-Xc)^2);
% 
% % https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwju_uvuhZXNAhVDE5QKHSQ8BqcQFggpMAI&url=http%3A%2F%2Fmy.ece.ucsb.edu%2FYork%2FBobsclass%2F194%2FReferences%2FGeneral%2FESR%2520of%2520Capacitors.pdf&usg=AFQjCNE8hpaUaAc0zW6mzUcYN4k5KOjT1g&bvm=bv.123664746,d.dGo&cad=rja
% ESR = Ras+1/(w^2*C^2*Rl)+1/(w^2*C^2*Rd)

%% Analyzing efficiency (power losses)(Q = around 100 in the kHz regime) for an inverted L-network
% Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf.
% See also "Insertion loss". pag 16 of lect7_match.pdf (https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwieq-Tao5fNAhXEW5QKHXKtAEgQFggdMAA&url=http%3A%2F%2Frfic.eecs.berkeley.edu%2F142%2Fpdf%2Flect7_match.pdf&usg=AFQjCNF-SM1zdrBem8kFllXPBTCAM1pQkg)
% Insertion loss = 1-efficiency;
f = 1e4;
Rs = 71.7;
L = 57.4e-3;
C1 = 4.8629e-09;
Qc1 = 100;
C2 = 4.5584e-08;
Qc2 = 100;
en_amp = 1.3e-9;
in_amp = 0.8e-12;
Rmatch = en_amp/in_amp;
Q = sqrt(Rmatch/Rs-1);
efficiency = (1-Q/Qc2)/(1+Q/Qc1); % eqs 16 & 18 from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf.
% for Q/Qc << 1 (high efficiency network) it can be simplified to
% efficiency2 = 1-Q/Qc2-Q/Qc1; % eqs 16 & 18 from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf. This agrees with lect7_match.pdf
G_ideal = sqrt(Rmatch/Rs); % from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf. eq.24
G_effective = G_ideal*sqrt(efficiency);
% Calculate the thevening impedance seen from the amplifier
Xl = 1j*2*pi*f*L;
Xc1 = 1/(1j*2*pi*f*C1);
Xc2 = 1/(1j*2*pi*f*C2);
Zth = (Rs+Xl+Xc1)*Xc2/(Rs+Xl+Xc1+Xc2);
%% Analyzing noise from the quality factor of the capacitors (Q = around 100 in the kHz regime) for an inverted L-network
% formulas for convertion from series to parallel
ESR1 = abs(Xc1)/Qc1;
% C2_p = C2*Qc2^2/(1+Qc2^2);
% Xc2_p = 1/(1j*2*pi*f*C2_p);
ESR2 = abs(Xc2)/Qc2; % <= % Q = abs(Xc)/ESR;
EPR2 = (1+Qc2^2)*ESR2;
% Noise introduced by non-ideal capacitors in an inverted L network
en_C1 = en_thermal(ESR1);
in_C2 = in_thermal(EPR2);
en_Rs = en_thermal(Rs);
en_total_in = sqrt((en_Rs*G_effective)^2+(in_amp*abs(Zth))^2+en_amp^2+(en_C1*G_effective)^2+(in_C2*abs(Zth))^2);

%% Analyzing noise For a T network:
f = 1e4;
Rs = 71.7;
L = 57.4e-3;
C1 = 4.8629e-09;
Qc1 = 100;
C2 = 4.5584e-08;
Qc2 = 100;
C3 = 1;
Qc3 = 100;
en_amp = 1.3e-9;
in_amp = 0.8e-12;
Rmatch = en_amp/in_amp;
% Designed from a given Q
Q = 5;
% for Q/Qc << 1 (high efficiency network) it can be simplified to
efficiency3 = 1-Q/Qc3-Q/Qc2-Q/Qc1; % eqs 16 & 18 from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf. This agrees with lect7_match.pdf
% G_ideal = sqrt(Rint/min(Rs,Rmatch));
G_ideal = sqrt(Q^2+1);
G_effective = G_ideal*sqrt(efficiency);
% Calculate the thevening impedances
Xl = 1j*2*pi*f*L;
Xc1 = 1/(1j*2*pi*f*C1);
Xc2 = 1/(1j*2*pi*f*C2);
Xc3 = 1/(1j*2*pi*f*C3);
Zth_T = Xc3+(Rs+Xl+Xc1)*Xc2/(Rs+Xl+Xc1+Xc2); % Thevening inpedance seen from the amplifier
Zth_L = (Rs+Xl+Xc1)*Xc2/(Rs+Xl+Xc1+Xc2); % Thevening inpedance seen from the C3 capacitor
%% Analyzing noise from the quality factor of the capacitors (Q = around 100 in the kHz regime) for an inverted L-network
% formulas for convertion from series to parallel
ESR1 = abs(Xc1)/Qc1;
% C2_p = C2*Qc2^2/(1+Qc2^2);
% Xc2_p = 1/(1j*2*pi*f*C2_p);
ESR2 = abs(Xc2)/Qc2; % <= % Q = abs(Xc)/ESR;
EPR2 = (1+Qc2^2)*ESR2;
ESR3 = abs(Xc3)/Qc3;
% Noise introduced by non-ideal capacitors in an inverted L network
en_C1 = en_thermal(ESR1);
in_C2 = in_thermal(EPR2);
en_C3 = en_thermal(ESR3);
en_Rs = en_thermal(Rs);
en_total_in = sqrt((en_Rs*G_effective)^2+(in_amp*abs(Zth_T))^2+en_amp^2+(en_C1*G_effective)^2+(in_C2*abs(Zth))^2+en_C3^2);