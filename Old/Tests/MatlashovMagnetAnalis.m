clear all, close all, clc
% Constants
T = 25 + 273.15; % Room temperature in Kelvin
kb = 1.3806488e-23; % J/K

%Dimentions of the coil
Din = 20e-3;
Dout = 90e-3;
W1 = 14e-3;
N =  1400;
W2 = (Dout -Din)/2;
D = ((Dout -Din)/2)+Din; % In meters!
d = 0.51e-3; % In meters!
Nmax = (pi/4)*(W1*W2)/(d^2);

R=20;
L=70e-3;
Ftrans = 70e3; %V/T

fc= R/(2*pi*L);
Area = (pi/4)*D^2; %*(pi/4);
Tib = Area*N/L;

% Rf = Ftrans/Tib;
Rf = logspace(2,9,8);
% Amplifier definition
en = 1.2e-9; % V/sqrt(Hz)
in = 0.8e-12; % A/sqrt(Hz)
% Output voltage per tesla or % Field transmision
B = 1; % One Tesla
Vout = B*Tib.*Rf;
% According to Tashiro 2006 'Optimizal..."
En = en; % Amplifier's voltage noise
R_Rf = (R*Rf)./(R+Rf);
Ei = in*R_Rf; % Current noise
Er =sqrt(4*kb*T*R_Rf); % Thermal noise
E = sqrt(En.^2+Ei.^2+Er.^2); %
Eout = E.*(1+Rf./R); %

% Sensitivity from Vout = B*Tib.*Rf;
Bmin = Eout./(Tib.*Rf);
% SN
SN = Vout./Eout;

% Determine the voltage created on the sides of the coil
w = 3300;
V = N*Area*w*B;

% Minimum noise if following paper and COMPUTING NOISE IN THE WRONG WAY!!!
BminTheo = E/Ftrans; % En dominates
BminTheo*1e15; % To convert it to femtotesla
% It gives 20ft/sqrt(Hz)!!! Like stated in the paper

Ftrans/V;