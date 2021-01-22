% DRAFT current vs voltage sensitivity modes
% From timofeeva et all, 2011, "Optimized mini search coil magnetometer suited to large BW applications"
clear all, close all, clc

Rb = 45; % Ohms
Lb = 47e-3; %Henrys
Cb = 60e-9; % Faradays
Seq = 1; % flux equivalent surface

% Thermal noise
ro = 1.7e-8;
T = 25 + 273.15; % 25 C in Kelvins
kb = 1.3806488e-23; % J/K
e_nRb = sqrt(4*kb*T*Rb); % Volts

%% Frequencies of interest
freq= logspace(0,9,40);
omega = 2*pi*freq;

%% Voltage analysis (equation derived from Eth=... (date 30/7/2015 in notebook))
dEthdBe = -1j*omega*Seq./(1-Lb*Cb*(omega.^2)+1j*omega*Rb*Cb);
sensE = e_nRb./dEthdBe;

%% Current analysis
dIndBe = -1j*omega*Seq./(Rb+1j*omega*Lb);
sensI = e_nRb./dIndBe;

% Visualize transfer functions and sensitivity limits
h=figure;
loglog(freq,abs(dEthdBe),'b'),hold on,
loglog(freq,abs(dIndBe),'r')
grid on, xlabel('Freq(Hz)'),ylabel('V/T or A/T'),title('Transfer function')
legend('Voltage','Current')
hold off
h=figure;
loglog(freq,abs(sensE),'b'),hold on,
loglog(freq,abs(sensI),'r')
grid on, xlabel('Freq(Hz)'),ylabel('T/Hz^1^/^2'),title('Noise floor (T/Hz^1^/^2)')
legend('Voltage','Current')
hold off