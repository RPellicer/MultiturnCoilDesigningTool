% GA optimization 
clear all, close all
nvars = 4;
%%
% --------INPUTS---------
%     N_v = x(1); % Loops per layer
%     N_l = x(2); % Layers
%     do = x(3);  % [METERS] wire outer diameter insulation included
%     di = x(4);  % [METERS] cooper outer diameter
% 'setGlobalVariables.m' is the source of fixed parameters
%% ######################################################################## 
% GLOBAL VARIABLES HAVE TO BE DEFINED IN 'setGlobalVariables.m'. The script 
% 'setGlobalVariables.m' allows setting the global variables within the
% workers of the parallel computation
%  ########################################################################
%% Linear inequality constraints of the form
% A[]*x[] <= b.
A = [0,0,-0.95,1]; %  -0.95*Dtot + 1*Dwire <= 0
b = 0;    %0.03 is the value of Rmax
IntCon = [1,2];
% Vector of lower bounds
lb = [1,1,0.0001,0.0001];
% Vector of upper bounds
ub = [50,50,0.006,0.003];
% % #########################################
% Start the parallel computing
if exist('poolobj')
    delete(poolobj)
end
poolobj = parpool;
% Define the options
options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',1e-17,'UseParallel', true); % The 'TolFun' should be set at least one order 
fitnessfcn = @Final_Tuned_solver;
% of magnitude smaller than the expected target value
% ,'InitialPopulation',[20,20,0.0004,0.0003]
% Call the optimization function
[x,fval,exitflag,output,population,scores] = ga(fitnessfcn,nvars,A,b,[],[],...
    lb,ub,@simple_nonlinear_constraint,IntCon,options);
% kill the parallel computing
delete(poolobj)

%% Show the results
setGlobalVariables(); % This function sets the global variables
N_v = x(1); % Loops per layer
N_l = x(2); % Layers
do = x(3);  % [METERS] wire outer diameter insulation included
di = x(4);  % [METERS] cooper outer diameter
Rin = Rout_glob - N_l*do;
dr = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
dz = 0.0; % azimut filling factor
R_i = Rin; % Coil internal radius    f0 = (max(freqs_glob)+min(freqs_glob))/2;
f0 = (freqs_glob(1)+ freqs_glob(end))/2;
bw = freqs_glob(end)-freqs_glob(1);
Visual_glob = 1; % plot results
[sensTot,sens,sensf0,freqs,C1,C2,Rs,Ls,Cs,emf_per_T,C1_ESR, C2_ESR] = sensTunedLinv(e_n_glob,i_n_glob,e_OutStagNoise_glob,N_v,N_l,R_i,di,do,dr,dz,f0,bw,Camp_glob,Ramp_glob,Qcap1ESR_glob,Qcap2ESR_glob,GainAmp_glob,Visual_glob);
disp(' ')
disp('#######################################################################################################################################################')
disp(['Estimated optimal coil has ' num2str(x(1)) ' loops per layer, ' num2str(x(2)) ' layers, with spacing of  ' num2str(x(3)) ' in between wires and cooper diameter of ' num2str(x(4)) ' (all in meters)']);
disp(['C1 = ' num2str(C1) ', C2 = ' num2str(C2)])
disp(['Peak sensitivity = ' num2str(min(sens)) ', Average sensitivity = ' num2str(sensTot)])
disp('#######################################################################################################################################################')
disp(' ')
% Data to use with Pspice simulator to compare results:
disp('Data to use into Pspice simulator to double check results:')
PspiceFreqDep = abs(mean(emf_per_T./(2*pi*freqs)));
disp(['Coil: Rs = ' num2str(Rs) ', Ls = ' num2str(Ls) ', Cs = ' num2str(Cs) ', PspiceFreqDep = ' num2str(PspiceFreqDep)])
disp(['Matching network: C1 = ' num2str(C1) ', C2 = ' num2str(C2), ', C1_ESR = ' num2str(C1_ESR) ', C2_ESR = ' num2str(C2_ESR)])
disp(' ')
