clear all, close all
%% L Matching network tomography (Rl>Rs?)
%
% ---Xs---------Xm----
%             |
%             |
%             Xt(Ct)
%             |
%             |            

% Define the properties of the coil
Ls_n = 57.4e-3;
Rs_n = 71.7;
Cs_n = 8.2e-12;
% Frequency of interest
f0_n = 1e4;
% Noise match the amplifier
en = 1.3e-9;
in = 0.8e-12;
% Calculate the matching apparent resistance
Rl_n = en/in;
% Rl_n = 10;
% Calculation of tuning capacitor
w0_n = 2*pi*f0_n;
Ct_simpleTuning = 1/(w0_n^2*Ls_n); % From Xc(w0) = Xl(w0) => (2j*pi*f0*Ct)^-1=-2j*pi*f0*Ls)
% Calculation of tuned impedance
syms Rs Ls Ct p w0 Rl real
syms Zt
% Derived from my calculations
% 1/Zt =1/Xc + 1/(Xl+Rs); => 1/Zt = 1j*w0*Ct+(1/(1j*w0*Ls+Rs))
% 1/Zt = (-(w0^2)*Ct*Ls+1j*w0*Ct*Rs+1)/(1j*w0*Ls+Rs)
% Zt = (1j*w0*Ls+Rs)/(-(w0^2)*Ct*Ls+1j*w0*Ct*Rs+1)
% Zt = (1j*w0*Ls+Rs)/((1-w0^2*Ct*Ls)+1j*w0*Ct*Rs)
% Removing imaginary numbers from the denominator to find the real and
% imaginary parts:
numer = expand((1j*w0*Ls+Rs)*((1-w0^2*Ct*Ls)-1j*w0*Ct*Rs));
denom = expand((1-w0^2*Ct*Ls)^2-(1j*w0*Ct*Rs)^2);
Zt = numer/denom;
Rt = real(Xt2);
Xt = imag(Xt2);
Ct_eq = nonzeros(solve(Rt==Rl,Ct));
Ct_n = eval(subs(Ct_eq, [w0,Rs,Ls,p,Rl],[w0_n,Rs_n,Ls_n,pi,Rl_n]));
% Calculate the total impedance seen from the matching capacitor
for k = 1:length(Ct_n)
    Zt_n(k,1) = eval(subs(Xt2, [w0,Rs,Ls,p,Ct],[w0_n,Rs_n,Ls_n,pi,Ct_n(k)]));
end
% Calculate the imaginary impedance that has to be compensated by the matching capacitor
for k = 1:length(Ct_n)
    Xt_n(k,1) = eval(subs(Xt, [w0,Rs,Ls,p,Ct],[w0_n,Rs_n,Ls_n,pi,Ct_n(k)]));
end
Cm = [];
for k = 1:length(Ct_n)
       Cm = [Cm; -1/(Xt_n(k)*w0_n)]; % Xc = -j/(wC)= 1/(jwC) => C = 1/(Xc*jw)
%        Lm2 = [Lm2; Zt2_n(k)/(w0_n)]; % Xl = jwL => L = Xl/jw
end
% Double check if they are capacitors
iscap_t = Ct_n>=0;
iscap_m = Cm>=0;
all_cap =(iscap_t+iscap_m)==2;
% Visualize the resistance against frequency
figure(1)
legendtext = [];
title('Plot of the impedances seen from the capacitor')
for k =1:length(all_cap)
    if (all_cap(k) == 1)
        Xt2_plot = eval(subs(Xt2, [Rs,Ls,p,Ct],[Rs_n,Ls_n,pi,Ct_n(k)]));
        f0_plotn=logspace(1,6,10000);
        w0_plotn=2*pi*f0_plotn;
        Xt2_plotn = eval(subs(Xt2_plot,w0,w0_plotn));
        semilogx(f0_plotn,real(Xt2_plotn))
        hold on
        grid on
        semilogx(f0_plotn,imag(Xt2_plotn))
        legendtext = [legendtext; ['Rt config nr.' num2str(k)]; ['Zt config nr.' num2str(k)]];
    else
        disp(['Solution nr.' num2str(k) ' not posible with only capacitors'])
    end
end
legend(legendtext)
hold off