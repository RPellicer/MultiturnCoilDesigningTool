%% Current magnetometer optimization
clear all, close all, clc
% Constants
T = 25 + 273.15; % Room temperature in Kelvin
kb = 1.3806488e-23; % J/K

%Dimentions of the coil
c = 0.2:0.2:3; % In CENTIMETERS!
d = 0.005:0.01:0.5; % In CENTIMETERS!
c_M = c'*ones(1,size(d,2));
d_M = ones(size(c,2),1)*d;
N = (pi/4)*(c_M.^2)./(d_M.^2);

% Inductance values
L_B = (25.49e-6*c_M.^5./d_M.^4)*1e-3;
L = L_B*(pi/4)^2;
% AC Resistance
R_B = 1.621e-5*c_M.^3./d_M.^4;
R = R_B;
% Cutoff frequency
fc = R./(2*pi*L);
% Ampere per tesla
Tib = 3.53e4*d_M.^2./c_M;

% Noise & amplification
% Amplifier definition
en = 0.9e-9; % V/sqrt(Hz)
in = 1e-12; % A/sqrt(Hz)
Rf = 1e6;
% Output voltage per tesla
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

for i = 1: size(E,1)
    for j = 1: size(E,2)
        if Ei(i,j)>En
            domin.source{i,j} = 'Ei';
            domin.value{i,j} = Ei(i,j);
        else
            domin.source{i,j} = 'En';
            domin.value{i,j} = En;
        end
        if Er(i,j)> domin.value{i,j}
            domin.source{i,j} = 'Er';
            domin.value{i,j} = Er(i,j);
        end
    end
end


% maximun wire thickness according to skin effect
f = [2000, 10000, 20000]; %hz
ro = 1.68e-8;
mu0 = pi*4e-7;
mur = 0.999994;

skinDepth = sqrt(2*ro./(2*pi*f*mu0*mur))*1e3