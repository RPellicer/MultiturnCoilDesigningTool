%% Sensitivity equation from Savukov 2007 test against data from Matlashov 2011
clc, clear cll, close all
% Bm = ;
% gamma = ;
freq = 3.3e3; % gamma*Bm; 
omega = 2*pi*freq;
d = 0.51e-3;
D_in = 20e-3;
D_out = 90e-3;
D = (D_out+D_in)/2;
W_1 = (D_out-D_in)/2;
W_2 = 14e-3;
N_t = 1400;
ro = 1.7e-8;
T = 25 + 273.15; % 25 C in Kelvins
kb = 1.3806488e-23; % J/K
ampNoise = 1.2e-9; % V/sqrt(Hz)

% Estimated number of turns N should be = N_t = 1400 /bad agreement 1400~4797/
N = 4*W_1*W_2/(pi*(d^2));

% field transfer coefficient (should give 70V/mT=7e4V/T)// good agreement 7e4=6.9e4 V/T  
vonB = omega*pi*(D^2)*N_t/4; % V/T
% v = B*omega*pi*(D^2)*N_t/4;

% Johnson noise & resistance
R_ac2 = N_t*D*pi*ro/(pi*((d/2)^2)); % Good agreement 20 = 20.1 ohms but it is DC resistance
V_rac2 =sqrt(4*kb*T*R_ac2); % 0.58nV/sqrt(Hz) estimated
% V_Rac3=sqrt(16*kb*T*ro*D*N_t / (d^2)); %

% sensitivity
V_w = pi*D*W_1*W_2; % coil 3D volume
coilSens = (8/(omega*D))*sqrt(kb*T*ro/V_w);% calculated 4.5fT/sqrt(Hz)~20fT/sqrt(Hz)reported
ampSens = (coilSens/V_rac2)*ampNoise;
totalSens = sqrt(coilSens^2+ampSens^2);% calculated 10fT/sqrt(Hz)~20fT/sqrt(Hz)reported

% MY SENSITIVITY EQUATION ESTIMATION
coilSensMy = V_rac2/vonB;
totalSensMy = sqrt(coilSensMy^2+ampSens^2);% calculated 10fT/sqrt(Hz)~20fT/sqrt(Hz)reported
