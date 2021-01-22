%% Calculating the sensitivity of the multiloop coil
% this function evaluates the sensitivity for known data points
% Number of turns per layer and Number of layers
% Diameter of wire with and without coil
% close all, clear all, clc

%Define global variables
global u0
u0=4*pi*1e-7; %permeability of free space
T = 25 + 273.15; % 25 C in Kelvins
k = 1.3806488e-23; % J/K

%% Calculate the teoretical sensitivity of the sensor
%% Electrical model
% Design by number of turns
di = 0.8e-3; % [METERS] cooper outer diameter
do = 0.83e-3; % [METERS] wire outer diameter insulation included

% To match indutance
N_v = 21; % Number of turns per layer
N_l = 5; % Number of layers
dr = 0.0e-3; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
dz = 0.0e-3; % azimut filling factor
s = do + sqrt(dr^2+dz^2); % distance between wire centres
R_i = (115/2)*1e-3; % Coil internal radius

% Noise from the amplifier (fixed to the amplifier type INA217)
i_n = 0.8e-12; %INA217
e_n = 1.3e-9;  %INA217
e_OutStagNoise = 90e-9;
Rf = 2*5000; % Feedback resistance of the amplifier (2x in the case if the INA)

% Frequencies of interest
freqs= 1e4;


sensTot = sensitivityI2V(i_n,e_n,e_OutStagNoise,Rf,N_v,N_l,R_i,di,do,dr,dz,freqs)
