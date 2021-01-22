clear all, close all, clc
%%  Design by number of turns
N_v = 35; % Number of turns per layer
N_l = 42; % Number of layers
do = 0.30e-3;      % [METERS] wire outer diameter insulation included
di = 0.5e-3;      % [METERS] cooper outer diameter
Rout = 18e-3; % Coil external radius;
Rin = Rout - N_l*do;
% To match indutance
dr = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
dz = 0.0; % azimut filling factor
%     s = do + sqrt(dr^2+dz^2); % distance between wire centres
R_i = Rin; % Coil internal radius

% Define the volume of the sample
r_ROI = 26e-3;       % radius of ROI
hBase_ROI = 4e-3;   % base height of ROI
hfinal_ROI = 14e-3;  % final height of ROI
Bp = 0.1; % Tesla 

%% Define constants
% Frequencies of interest
% freqs = 3.3e3;
freqs = logspace(2,5,20);
% Noise from the amplifier (fixed to the amplifier type INA217)
i_n = 0.8e-12;          %INA217
e_n = 1.2e-9;           %INA217
e_OutStagNoise = 90e-9; %INA217
%   % Ideal amplifier
%         i_n = 0.8e-20; %INA217
%         e_n = 1.2e-20;  %INA217
%         e_OutStagNoise = 90e-20;
% Feedback resistance of the amplifier (2x in the case if the INA)
Rf = 2*5000; 

[stotal,sens] = sensitivityI2V(i_n,e_n,e_OutStagNoise,Rf,N_v,N_l,R_i,di,do,dr,dz,freqs);
[SNR_ideal,vn_in,emf] = SNR_volume(i_n,e_n,e_OutStagNoise,Rf,N_v,N_l,R_i,di,do,dr,dz,freqs,r_ROI,hBase_ROI,hfinal_ROI,Bp);
% Account for prepolarization efficiency, AC character and exponential
% decay when calculating the expected emf in the espectra
BpTime = 10; % Ten seconds of prepolarization
relaxTimeBp = 3; % T1 = 3 seconds for the water @ 0.1 T
relaxTimeBm = 1; % T1 = 12 seconds for the water @ 50 uT
t_0 = 20e-3; 
t_end = 1;
emf_prep = emf*(1-exp(-BpTime/relaxTimeBp));
emf_eff = MeanExpEnvelope(emf_prep,relaxTimeBm,t_0,t_end)/sqrt(2);
% Calculate the realistic SNR
SNR = abs(emf_eff./vn_in);

% Visualize
% The sensitivity
figure(1)
loglog(freqs,sens,'b')
hold on
line([freqs(1),freqs(end)],[stotal,stotal],'Color','r')
grid on, title('Sensitivity (T/Hz^1^/^2)'), legend('Sensitivity', 'Mean Sensitivity')
hold off
% The SNR
figure(2)
loglog(freqs,SNR_ideal,'b')
hold on,
loglog(freqs,SNR,'r')
grid on, title('SNR'), hold off
legend('Ideal DC SNR', 'AC,exp and prepolarization efficiency accounted SNR')
% The Noise
figure(3)
loglog(freqs,abs(vn_in),'b')
hold on,
% The EMF
loglog(freqs,abs(emf),'r')
loglog(freqs,abs(emf_prep),'g')
loglog(freqs,abs(emf_eff),'k')
legend('noise','ideal emf_p_p','Prepolarization time accounted emf','Prep time, sinusoidal and exponential decay accounted emf')
grid on, title('Noise & Signal (volts)'), hold off
