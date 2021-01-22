%% GA optimization 
% ######################################################################## 
% GLOBAL VARIABLES HAVE TO BE DEFINED IN 'setGlobalVariables.m'. The script 
% 'setGlobalVariables.m' allows setting the global variables within the
% workers of the parallel computation
% ########################################################################
clear, close all
setGlobalVariables()
% setGlobalVariablesHerrer()

% ########################################################################
% Sensing method:  Non-tuned    Tuned        Both
%     type =          0           1           2
type = 0;       % ######## ADJUSTABLE ###########
maxCount = 4;   % ######## ADJUSTABLE ###########
useparallel = 1;% ######## ADJUSTABLE ###########
% Visualize
graphs = 0; % plot results
seeProfile = 0;
% Timing
procesHours = 20;
procesMins = 5;
% ########################################################################

%% Initialize
% run more that once the optimization and archive the results
disp(['Computation started: ' datestr(now)]);
tic
elapsedTime = toc;
count = 1;
sweepCount = 1;
%% Set the properties of amplifier, coil and matching network elements
% --------INPUTS---------
%     N_v = x(1); % Loops per layer
%     N_l = x(2); % Layers
%     do = x(3);  % [METERS] wire outer diameter insulation included
%     di = x(4);  % [METERS] cooper outer diameter
nvars = 2;
exitflag = 10; % The 'ga' will overwrite this value to -1 when the optimization is stopped by the user. The value of 10 has been given because values from -5 to 5 have their own meaning in 'ga'
%% Linear inequality constraints of the form
% A[]*x[] <= b.
A = [-0.95 1]; %  -0.95*Dtot + 1*Dwire <= 0
b = 0;    %0.03 is the value of Rmax
% IntCon = [1,2];
% Vector of lower bounds
lb = [0.0001 0.0001];
% Vector of upper bounds
ub = [0.005 0.005*0.95];
% Lower bound on the change in the value of the objective function during a step
if tProbe.Coil.emfNear
    tolerFun = 1e-5;
else
    tolerFun = 5e-16;
end
% #########################################
% Start the parallel computing
if useparallel
    if exist('poolobj')
        delete(poolobj);
    end
    poolobj = parpool;
end
%% Prepare counters and flags for the loops
cType = type; % Here it sets tuned to 0, 1 or 2
countType = 1;
maxCountType = 1;
if type == 2
    cType = 1; % If type ==2, sets tuned to 0 first and then it alternates its value
    maxCountType = 2;
end
ProbeOpt = optimProbe; % "optimProbe" is a class
%% Loop
while (toc < (procesHours*60^2 + procesMins*60)) && (count<=maxCount) && (sweepCount<=maxSweepCount)  && (exitflag ~= -1)
%     tProbe.Coil.LitzStrandNr_glob = sweepOptions(sweepCount); % Sweep this value for optimization process
    tProbe.(tProbe.Freqs.sweepPropA).(tProbe.Freqs.sweepPropB) = sweepOptions(sweepCount); % Sweep this value for optimization process
%     tProbe.MatchNetw.Qcap2ESR_glob = sweepOptions(sweepCount);
    tProbe.Coil.Lmax_glob = sweepOptions(sweepCount)/2;  % The maximum lift of a Brooks coil
    %% Optimization process
    if cType == 1
    % Define the options: The 'TolFun' should be set at least one order of magnitude smaller than the expected target value
    % Optimization terminated: average change in the penalty fitness value less than options.TolFun and constraint violation is less than options.TolCon.
        options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',tolerFun,'UseParallel', useparallel); 
        fitnessfcn = @(x)Final_Tuned_solver_Brooks(x,tProbe,Visual_glob);
        tProbe.Type = 'Tuned';
    else
        options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',tolerFun,'UseParallel', useparallel);
        fitnessfcn = @(x)Final_solver_Brooks(x,tProbe,Visual_glob);
        tProbe.Type = 'Non-tuned';
    end
    % ,'InitialPopulation',[20,20,0.0004,0.0003]
    % Call the optimization function
    if seeProfile
        profile on
    end
    [x,fval,exitflag,output,population,scores] = ga(fitnessfcn,nvars,A,b,[],[],...
        lb,ub,[],options);
    %% Show the results
    if seeProfile
        profile off
        profile viewer
    end
    tProbe.Coil.do = x(1);  % [METERS] wire outer diameter insulation included
    tProbe.Coil.di = x(2);  % [METERS] cooper outer diameter
    tProbe.Coil.N_v = floor(tProbe.Coil.Rout_glob/(2*tProbe.Coil.do));
    tProbe.Coil.N_l = tProbe.Coil.N_v;
    % Set other needed variables
    tProbe.Coil.Rin = tProbe.Coil.Rout_glob - tProbe.Coil.N_l*tProbe.Coil.do;
    tProbe.Sens.OptProcTime = toc - elapsedTime;
    % #####################################################################################################################################################################
    disp(['Iteration Nr.: ' num2str(countType + ((sweepCount-1)*maxCountType) + ((count-1)*maxSweepCount*maxCountType))])  
    disp('#######################################################################################################################################################')
    StrProbeOpt = getProbeDetails(tProbe,graphs);
    ProbeOpt(count,sweepCount,countType).Type = StrProbeOpt.Type;
    ProbeOpt(count,sweepCount,countType).Coil= StrProbeOpt.Coil;
    ProbeOpt(count,sweepCount,countType).Amp= StrProbeOpt.Amp;
    ProbeOpt(count,sweepCount,countType).Sens = StrProbeOpt.Sens;
    ProbeOpt(count,sweepCount,countType).Freqs = StrProbeOpt.Freqs;
    if cType == 1
        ProbeOpt(count,sweepCount,countType).MatchNetw = StrProbeOpt.MatchNetw;
    end
    % ######################################################################################################################################################################
    elapsedTime = toc;
    save('optimResults.mat','ProbeOpt','elapsedTime'); % Saves the last data in case it crashes
    % Update counters
    countType = countType +1;
    if countType > maxCountType
        countType = 1;
        sweepCount = sweepCount +1;
        if sweepCount > maxSweepCount
            sweepCount = 1;
            count = count +1;
        end
    end
    if type == 2 && countType ==1
        cType = 1; % If type ==2, sets tuned to 0 first and then it alternates its value
    elseif type == 2 && countType ==2
        cType = 0; % If type ==2, sets tuned to 0 first and then it alternates its value
    end
    % #####################################################################################################################################################################
end
% Stop the parallel computing
if exist('poolobj')
    delete(poolobj);
end
% Save the data
save('optimResults.mat','ProbeOpt','elapsedTime')
% This will shut down the computer in 3 minutes time. It can be reversed by >>
% system('shutdown -a')
system('shutdown -s -t 180')
% VisualizeCoilOptRes()
% To visualize one of the results:
% getProbeDetails(ProbeOpt(3,1,1),1);