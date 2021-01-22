%% #########################################
% Global variables to be adjusted from here. All values are International
% System of Units unless otherwise specified.

%% For sweeping. IMPORTANT: The values defined in "sweepOptions" will OVERWRITE the values of the parameter.
% First we define the parameter that needs to be sweept
tProbe.Freqs.sweepPropA = 'Freqs';  % Change here the variable that you want to sweep.
% tProbe.Freqs.sweepPropA = 'Coil';  
% tProbe.Freqs.sweepPropA = 'MatchNetw';  
% tProbe.Freqs.sweepPropB = 'f0_glob';
% tProbe.Freqs.sweepPropB = 'bw_glob'; % Change here the variable that you want to sweep.
tProbe.Freqs.sweepPropB = 'f0_Offset'; % Change here the variable that you want to sweep.
% tProbe.Freqs.sweepPropB = 'LitzStrandNr_glob';
% tProbe.Freqs.sweepPropB = 'Rout_glob';
% Then we define the values to be used for that parameter
% sweepOptions = [-200 -100 0 100 200]; % Set to '0' if no parameter needs to be sweept
sweepOptions = [0];
% sweepOptions = logspace(1,log10(19900),10);
% sweepOptions = [3.3e3 10e3]; % f0
% sweepOptions = logspace(log10(1e3),log10(100e3),5); % f0
% sweepOptions = [9700 9800 9900 10000 10100 10200 10300] - 10e3; % f0_Offset
% sweepOptions = [3100 3300 3500] - 3.3e3; % f0_Offset
% sweepOptions = [33 100 330 1000 1000000000]; % 33 100 330 1000]; Q factor 
% sweepOptions = [0.01 0.015 0.02 0.025 0.03]; % Rout for BROOKs 0.01 0.02 0.03 0.04 0.05
% litzOptions = round(logspace(0,3,8));
% tProbe.Freqs.f0_glob = logspace(log10(1e3),log10(500e3),20);% 100e3;
tProbe.Coil.LitzStrandNr_glob = 1; % Set to 1 for normal wire (no Litzwire)
% Coil Shape    
% tProbe.Coil.Rout_glob = 90e-3;  % Outer radious of the coil (fixed)
% tProbe.Coil.Lmax_glob = 18e-3;  % Maximum height of the coil
%% 50ml Sample holder
% tProbe.Coil.Rout_glob = 45e-3;  % Outer radious of the coil (fixed)
% tProbe.Coil.Lmax_glob = 40e-3;  % Maximum height of the coil % 20ml in the sample holder takes 28.5e-3 length
tProbe.Coil.Rout_glob = (54e-3)/2 + 0.8e-3;  % Outer radious of the coil (fixed)
tProbe.Coil.Lmax_glob = 80e-3;  % Maximum height of the coil % 20ml in the sample holder takes 28.5e-3 length
% These are not used so are set to 0
tProbe.Coil.dr = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
tProbe.Coil.dz = 0.0; % azimut filling factor
% Others
tProbe.Coil.emfNear = 1; % If sensitivity to a specific volume in the near field is needed set this to 1. CHECK NEXT 4 PARAM
tProbe.Coil.Rin_glob = (54e-3)/2; % 45/2; %You can define the inner diameter of the coil for "tProbe.Coil.emfNear".
tProbe.Coil.Bp = 45e-3;                 % Prepolarization field strength (This helps calculating the SNR)/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.r_ROI = (50e-3)/2;% tProbe.Coil.Rout_glob/2;    % RADIUS of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.hBase_ROI = -25e-3;         % Closes part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.hfinal_ROI = 25e-3;%tProbe.Coil.Rout_glob; % Furthest part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.Coil_heigth_cent = 1; % Only for emfNear = 1 % tProbe.Coil.Coil_heigth_centered at location 0 mm = 1;
tProbe.Coil.PenalyzeScaling = 1;
% Matching network
tProbe.MatchNetw.Qcap1ESR_glob = 1000;
tProbe.MatchNetw.Qcap2ESR_glob = 1000;
% Frequencies of interest
tProbe.Freqs.bw_glob = 5e3; % Bandwidth of interest
tProbe.Freqs.f0_glob = 78e3; % Central frequency of interest
tProbe.Freqs.f0_Offset = 0;
% Amplifier
% #############################################
% % INA217
% tProbe.Amp.Camp_glob = 2e-12;
% tProbe.Amp.Ramp_glob = 60e6;
% tProbe.Amp.GainAmp_glob = 1000;
% tProbe.Amp.Rf_glob = 2*5000;
% % NOISE
NrAmpParall = 1; % Number of amplifiers in parallel. Reduces voltage noise but increases current noise
% tProbe.Amp.e_n_glob = 1.3e-9/sqrt(NrAmpParall);
% tProbe.Amp.i_n_glob = 800e-15/sqrt(NrAmpParall);
% tProbe.Amp.e_OutStagNoise_glob = 90e-9*sqrt(NrAmpParall); % en_thermal(tProbe.Amp.Rf_glob);
% % LNA718
% tProbe.Amp.e_n_glob = 0.7e-9/sqrt(NrAmpParall);
% tProbe.Amp.i_n_glob = 98e-15/sqrt(NrAmpParall);
% tProbe.Amp.e_OutStagNoise_glob = en_thermal(tProbe.Amp.Rf_glob)*sqrt(NrAmpParall);
% % % My test
% tProbe.Amp.e_n_glob = 0.5e-9/sqrt(NrAmpParall);
% tProbe.Amp.i_n_glob = 80e-15/sqrt(NrAmpParall);
% tProbe.Amp.e_OutStagNoise_glob = 1e-30;
% % LNAMFBX-1000
% tProbe.Amp.Camp_glob = 56e-12;
% tProbe.Amp.Ramp_glob = 100e3;
% tProbe.Amp.GainAmp_glob = 1000;
% tProbe.Amp.Rf_glob = 2*5000;
% tProbe.Amp.e_n_glob = 0.53e-9/sqrt(NrAmpParall);
% tProbe.Amp.i_n_glob = 15e-15/sqrt(NrAmpParall);
% tProbe.Amp.e_OutStagNoise_glob = 1e-30;
% % Miguel-Janasek 8x 2SK3557
tProbe.Amp.Camp_glob = 56e-12;
tProbe.Amp.Ramp_glob = 100e6;
tProbe.Amp.GainAmp_glob = 1000;
tProbe.Amp.Rf_glob = 2*500000;
tProbe.Amp.e_n_glob = 1.05e-9/sqrt(NrAmpParall); % 0.44e-9/sqrt(NrAmpParall);
tProbe.Amp.i_n_glob = (4.4e-15)*sqrt(NrAmpParall);
tProbe.Amp.e_OutStagNoise_glob = 1e-30;
% % Ideal
% tProbe.Amp.e_n_glob = 1e-30;
% tProbe.Amp.i_n_glob = 1e-30;
% tProbe.Amp.e_OutStagNoise_glob = 1e-30;
% LT1028 Tashiro-2006-Opt...
% tProbe.Amp.e_n_glob = 0.9e-9/sqrt(NrAmpParall);
% tProbe.Amp.i_n_glob = 100e-15/sqrt(NrAmpParall);
% tProbe.Amp.e_OutStagNoise_glob = en_thermal(tProbe.Amp.Rf_glob)*sqrt(NrAmpParall);

Visual_glob = 0; % Visualize every calculation? plot results (1) or not (0)
maxSweepCount = length(sweepOptions);
tProbe.Coil.discreteWire = 1;