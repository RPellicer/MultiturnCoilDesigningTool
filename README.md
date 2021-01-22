The main file you should execute is “Main_GA_Solver_flexible.m “.

Please, adjust the parameters as needed. In the code they are set to:

·         Coil outer radius = 25e-3.

·         Maximum coil height=10e-3.

·         Quality factor of the capacitors = 200. You can set the quality factor of the capacitors manually. Although for the paper I chose one capacitor series and characterised it, I personally find it more practical to set this value by hand to the type of capacitors you are using.

·         Central frequency = 3e3.

·         Frequency bandwidth = 100.

·         Preamplifier INA217.

These parameters are defined in “setGlobalVariables.m” such that:

…

% Coil Shape   
tProbe.Coil.Rout_glob = 25e-3;  % Outer radious of the coil (fixed)
tProbe.Coil.Lmax_glob = 10e-3;  % Maximum height of the coil

…

% Matching network

tProbe.MatchNetw.Qcap1ESR_glob = 200;

tProbe.MatchNetw.Qcap2ESR_glob = 200;

% Frequencies of interest

tProbe.Freqs.bw_glob = 100; % Bandwidth of interest

tProbe.Freqs.f0_glob = 3e3; % Central frequency of interest

% Amplifier

% #############################################

% INA217

…

% NOISE

…

tProbe.Amp.e_n_glob = 1.3e-9/sqrt(NrAmpParall);

tProbe.Amp.i_n_glob = 800e-15/sqrt(NrAmpParall);

tProbe.Amp.e_OutStagNoise_glob = 90e-9*sqrt(NrAmpParall); % en_thermal(tProbe.Amp.Rf_glob);

…

 

Other things you may want to play with:

·         Check the sensitivity for non-tuned current-to-voltage vs. tuned voltage-to-voltage amplification by setting in “Main_GA_Solver_flexible.m “:

type = 2;       % ######## ADJUSTABLE ###########

·         Run the optimisation a number of times by changing in “Main_GA_Solver_flexible.m “:

maxCount = 3;

·         Speed up the process by parallel using multiple processors (parallel CPUs).

useparallel = 1;

·         Sweep one parameter to see how it affects the sensitivity (e.g., frequency bandwidth or coil outer diameter) by setting in “setGlobalVariables.m” the name of the parameter. Because of the architecture of the data I have to define 2 names: For example, to sweep the outer radius of your coil.

tProbe.Freqs.sweepPropA = 'Coil';  % Change here the variable that you want to sweep.

tProbe.Freqs.sweepPropB = 'Rout_glob'; % Change here the variable that you want to sweep.

sweepOptions = [0.01 0.015 0.02 0.025 0.03];

 

An important parameter that may need to be adjusted is “tolerFun” in “Main_GA_Solver_flexible.m”. I usually set it to one or two orders of magnitude smaller than the highest sensitivity I would expect from the detector. This basically tells the program: when your search doesn’t improve the sensitivity by more than that value for a while, stop searching. I would set it to 1e-15 in your case. While the optimisation algorithm is running you will see a doted plot being updated. The value you see under “Best” is the highest mean sensitivity it has a found. If you set the “tolerFun” too demanding, i.e., too small, the optimisation program may keep trying to search for a long time. You can stop it anytime by pressing on the “stop” button shown on the lower left side of the plot window. The computer may need some time to react to this if you are using all your CPUs (parallel computing option).

 

If you find the optimal values given by the program out of your available range, such as thinner wire than the one you can find, you will need to adjust the boundary conditions by changing in “Main_GA_Solver_flexible.m“:

% Vector of lower bounds

lb = [1 1 0.000101 0.0001]; % minimum [Number of loops per layer (N_v), Number of layers (N_l), Wire outer diameter (do), Wire cooper outer diameter (di)]

% Vector of upper bounds

ub = [80 80 0.005 0.005*0.95]; % maximum [Number of loops per layer (N_v), Number of layers (N_l), Wire outer diameter (do), Wire cooper outer diameter (di)]

 

I find re-running the optimisation several times very useful (e.g., maxCount = 6;). As you may observe, there are several configurations that reach similar sensitivities. Some of them will be more practical to build. In my case, I pick the one that most closely matches commercial wires.

 

The results from the optimisation will be displayed on Matlab’s command window as they are being calculated. This results are stored under the object ProbeOpt. If you want to visualise one specific result you can by the function “getProbeDetails(ProbeOpt(index),1)”. You will need to give the index that corresponds to the coil you want to visualise, e.g., getProbeDetails(ProbeOpt(1,1,1),1). This will visualise a set of useful plots and information, such as the noise contribution of each component (as shown in the paper).

 

Another option is to test certain configuration manually, which can be done by “CalcSensFromDesign(NrLoopsPerLayer,NrLayers,OutWireDiam,CoopDiam, AmplMode)”. This is useful because it allows us testing configurations with off-the-shelf components in order to see if the sensitivity lost is significant as compared with the ideal ones.

 

Hope that the program helps you achieving the sensitivity you need for your application.

Don’t hesitate getting back to me if you have any question.
