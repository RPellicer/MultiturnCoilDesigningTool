# Multiturn coil designing tool
A tool to design low frequency (> ~100 kHz) air-core magnetometers

## How to execute it
The code is run by executing "Main_GA_Solver_flexible.m" from MatLab.

Please, adjust the parameters as needed. These parameters are defined in "setGlobalVariables.m" such that:

```MatLab
%% #############   General   ############################
tProbe.Coil.discreteWire = 0; % 0 = any value, 1 = integer values (steps of 0.1 mm for wire and coating sizes to accelerate the search)

%% #############   sweeping a parameter   ############################
% For sweeping a parameter. IMPORTANT: The values defined in "sweepOptions" will OVERWRITE the values of the parameter.
% First we define the parameter that needs to be sweept
tProbe.Freqs.sweepPropA = 'Freqs';  % Change here the variable that you want to sweep.
tProbe.Freqs.sweepPropB = 'f0_Offset'; % Change here the variable that you want to sweep.
sweepOptions = [0]; % Set this to '0' if you don't want to sweep any parameter
% sweepOptions = [-200 -100 0 100 200]; % Set to '0' if no parameter needs to be sweept

%% ###############  Coil  ##########################
%% Litz-wire
tProbe.Coil.LitzStrandNr_glob = 1; % Set to 1 for normal wire (not Litzwire)
tProbe.Coil.NS = tProbe.Coil.LitzStrandNr_glob; % total number of strands. Only relevant if tProbe.Coil.LitzStrandNr_glob > 1
tProbe.Coil.NB = 1;                 % number of bunching operations, (agrupations of strands). Only relevant if tProbe.Coil.LitzStrandNr_glob >
tProbe.Coil.NC = 1;                 % number of cabling,(agrupations of bunches). Only relevant if tProbe.Coil.LitzStrandNr_glob >

% Sensitivity to a volume in the near field
tProbe.Coil.emfNear = 0; % If sensitivity to a specific volume in the near field is needed set this to 1. CHECK NEXT 4 PARAM
tProbe.Coil.Bp = 50e-3;                 % Prepolarization field strength (This helps calculating the SNR)/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.r_ROI = (50e-3)/2;% tProbe.Coil.Rout_glob/2;    % RADIUS of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.hBase_ROI = -25e-3;         % Closes part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.hfinal_ROI = 25e-3;%tProbe.Coil.Rout_glob; % Furthest part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
tProbe.Coil.Coil_heigth_cent = 1; % Only for emfNear = 1 % tProbe.Coil.Coil_heigth_centered at location 0 mm = 1;

% Coil Shape contraints    
tProbe.Coil.Lmax_glob = 20e-3;  % Maximum height of the coil
tProbe.Coil.Rout_glob = 45e-3;  % Outer radious of the coil (fixed), not used if tProbe.Coil.emfNear = 1
tProbe.Coil.Rin_glob = 45/2; % not used if tProbe.Coil.emfNear = 0

% These are not used so are set to 0
tProbe.Coil.dr = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
tProbe.Coil.dz = 0.0; % azimut filling factor

% General
tProbe.Coil.PenalyzeScaling = 1;

%% ###############  Coil  ##########################
% Matching network
tProbe.MatchNetw.Qcap1ESR_glob = 1000;
tProbe.MatchNetw.Qcap2ESR_glob = 1000;
% Frequencies of interest
tProbe.Freqs.bw_glob = 1e3; % Bandwidth of interest
tProbe.Freqs.f0_glob = 10e3; % Central frequency of interest
tProbe.Freqs.f0_Offset = 0;

%% Amplifier
% #############################################
% % INA217
tProbe.Amp.Camp_glob = 2e-12;
tProbe.Amp.Ramp_glob = 60e6;
tProbe.Amp.GainAmp_glob = 1000;
tProbe.Amp.Rf_glob = 2*5000;
% NOISE
NrAmpParall = 1; % Number of amplifiers in parallel. Reduces voltage noise but increases current noise
tProbe.Amp.e_n_glob = 1.3e-9/sqrt(NrAmpParall);
tProbe.Amp.i_n_glob = 800e-15/sqrt(NrAmpParall);
tProbe.Amp.e_OutStagNoise_glob = 90e-9*sqrt(NrAmpParall); % en_thermal(tProbe.Amp.Rf_glob);

%% Visualise everything (not recomended)
Visual_glob = 0; % Visualize every calculation? better set "graph" to 1 on "Main_GA_Solver_flexible.m"
maxSweepCount = length(sweepOptions);
```

## Other adjustments

An important parameter that may need to be adjusted is "tolerFun" in "Main_GA_Solver_flexible.m". I usually set it to one or two orders of magnitude smaller than the highest sensitivity I would expect from the detector. This basically tells the program: when your search doesn't improve the sensitivity by more than that value for a while, stop searching. You can run the program once to get an estimate of what sensitivity you can expect. While the optimisation algorithm is running you will see a doted plot being updated. The value you see under "Best" is the highest mean sensitivity it has a found. If you set the "tolerFun" too demanding, i.e., too small, the optimisation program may keep trying to search for a long time. You can stop it anytime by pressing on the "stop" button shown on the lower left side of the plot window. The computer may need some time to react to this if you are using all your CPUs (parallel computing option).

If you find the optimal values given by the program out of your available range, such as thinner wire than the one you can find, you will need to adjust the boundary conditions by changing in "Main_GA_Solver_flexible.m":
```MatLab
% Vector of lower bounds
lb = [1 1 0.00001 0.0001]; % minimum [Number of loops per layer (N_v), Number of layers (N_l), free space betwen wires (do-di)(the couting can be used some times to materialise this spacing), Wire cooper outer diameter (di)]
% Vector of upper bounds
ub = [140 15 0.015 0.005*0.95]; % maximum [Number of loops per layer (N_v), Number of layers (N_l), free space betwen wires (do-di)(the couting can be used some times to materialise this spacing), Wire cooper outer diameter (di)] Lower bound on the change in the value of the objective function during a step
```
 
I find re-running the optimisation several times very useful (e.g., maxCount = 6;). As you may observe, there are several configurations that reach similar sensitivities. Some of them will be more practical to build. In my case, I pick the one that most closely matches commercial wires.

The results from the optimisation will be displayed on Matlab's command window as they are being calculated. This results are stored under the object ProbeOpt. If you want to visualise one specific result you can by the function "getProbeDetails(ProbeOpt(index),1)". You will need to give the index that corresponds to the coil you want to visualise, e.g., "getProbeDetails(ProbeOpt(1,1,1),1)". This will visualise a set of useful plots and information, such as the noise contribution of each component (as shown in the paper).

Another option is to test certain configuration manually, which can be done by "CalcSensFromDesign(NrLoopsPerLayer,NrLayers,OutWireDiam,CoopDiam, AmplMode)". This is useful because it allows us testing configurations with off-the-shelf components in order to see if the sensitivity lost is significant as compared with the ideal ones.

Hope that the program helps you achieving the sensitivity you need for your application.

Don't hesitate getting back to me if you have any question.

## For attribution, please cite this work as:
Pellicer-Guridi, R., Vogel, M.W., Reutens, D.C. et al. Towards ultimate low frequency air-core magnetometer sensitivity. Sci Rep 7, 2269 (2017). https://doi.org/10.1038/s41598-017-02099-z
