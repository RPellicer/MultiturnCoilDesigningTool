% SNR estimation without losses
clear all, close all
Rmatch = logspace(-1,7,1000);
Rs = 0.7;
en = 1.3e-9;
in = 0.8e-12;
m = Rmatch./Rs;% Rmatch./Rs=(Vin/Vout)^2
Gain = sqrt(m); % The voltage gain is defined by the resistance convertion 'm'
Zth = Rmatch;
SNR = Gain./sqrt(en^2+(in.*Zth).^2+(en_thermal(Rs)*Gain).^2);
semilogx(Rmatch,SNR)
grid on