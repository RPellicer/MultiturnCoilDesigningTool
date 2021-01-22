%% GA optimization 
% ######################################################################## 
% GLOBAL VARIABLES HAVE TO BE DEFINED IN 'setGlobalVariables.m'. The script 
% 'setGlobalVariables.m' allows setting the global variables within the
% workers of the parallel computation
% ########################################################################
clear, %close all
setGlobalVariables()
% setGlobalVariablesHerrer()

% ########################################################################
% Sensing method:  Non-tuned    Tuned        Both
%     type =          0           1           2
type = 2;       % ######## ADJUSTABLE ###########
maxCount = 2;   % ######## ADJUSTABLE ########### Rerun the same optimisation "maxCount" times to double check reproducivility of results.
useparallel = 1;% ######## ADJUSTABLE ########### Set to '1' for parallel processing. This may slow down other tasks you are doing with the computer
% Visualize
graphs = 1;     % plot results
seeProfile = 0; % insight into bottle neck functions which could be optimised for speed
% Time limiter for the whole process
turn_off_PC = 0;    % set to 1 to shut-down the PC after completion. Results in optimProbe
procesHours = 200;
procesMins = 5;
% ########################################################################
%% Initialize
% run more than once the optimization and archive the results
disp(['Computation started: ' datestr(now)]);
tic
elapsedTime = toc;
count = 1;
sweepCount = 1;
%% Set the properties of amplifier, coil and matching network elements
% --------INPUTS---------
%     N_v = x(1); % Loops per layer
%     N_l = x(2); % Layers
%     if tProbe.Coil.discreteWire == 1 % Increase wire and coating sizes in steps of 0.1 mm to accelerate the search
%     do = (x(4)+x(3))*1e-4;  % [METERS] wire outer diameter insulation included
%     di = x(4)*1e-4;  % [METERS] cooper outer diameter
%     if tProbe.Coil.discreteWire == 0;
%     do = x(3);  % [METERS] wire outer diameter insulation included
%     di = x(4);  % [METERS] cooper outer diameter
nvars = 4;
exitflag = 10; % The 'ga' will overwrite this value to -1 when the optimization is stopped by the user. The value of 10 has been given because values from -5 to 5 have their own meaning in 'ga'
%% Linear inequality constraints of the form A[]*x[] <= b.
if tProbe.Coil.discreteWire == 1
    A = [0 0 -100 1]; %  -100*do + 1*di <= 0 % To ensure a minimum insolation coating
    b = 0;    %0.03 is the value of Rmax    
    IntCon = [1,2,3,4]; % Which variables are integer [N_v, N_l], Ruben2020 to search quicker    
    lb = [1 1 1 1]; % minimum [Number of loops per layer (N_v), Number of layers (N_l), Wire outer diameter (do), Wire cooper outer diameter (di)] Ruben2020 to search quicker
    ub = [140 15 50 50]; % maximum [Number of loops per layer (N_v), Number of layers (N_l), Wire outer diameter (do), Wire cooper outer diameter (di)] Ruben2020 to search quicker
else
    A = [0 0 -0.95 1]; %  -0.95*do + 1*di <= 0 % To ensure a minimum insolation coating Ruben 2020
    b = 0;    %0.03 is the value of Rmax
    IntCon = [1,2]; % Which variables are integer [N_v, N_l]
    % Vector of lower bounds
    lb = [1 1 0.000101 0.0001]; % minimum [Number of loops per layer (N_v), Number of layers (N_l), Wire outer diameter (do), Wire cooper outer diameter (di)]
    % Vector of upper bounds
    ub = [140 15 0.005 0.005*0.95]; % maximum [Number of loops per layer (N_v), Number of layers (N_l), Wire outer diameter (do), Wire cooper outer diameter (di)] Lower bound on the change in the value of the objective function during a step
end
if tProbe.Coil.emfNear
    tolerFun = 1e-15;
else
    tolerFun = 1e-17;
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
while (toc < (procesHours*60^2 + procesMins*60)) && (count<=maxCount) && (sweepCount<=maxSweepCount) && (exitflag ~= -1)
%     tProbe.Coil.LitzStrandNr_glob = sweepOptions(sweepCount); % Sweep this value for optimization process
    tProbe.(tProbe.Freqs.sweepPropA).(tProbe.Freqs.sweepPropB) = sweepOptions(sweepCount); % Sweep this value for optimization process
%     tProbe.MatchNetw.Qcap2ESR_glob = sweepOptions(sweepCount);
%     tProbe.Coil.Lmax_glob = sweepOptions(sweepCount)/2;  % The maximum lift of a Brooks coil
%     tProbe.Freqs.bw_glob = 2*sweepOptions(sweepCount)*0.95;
    %% Optimization process
    if cType == 1
    % Define the options: The 'TolFun' should be set at least one order of magnitude smaller than the expected target value
    % Optimization terminated: average change in the penalty fitness value less than options.TolFun and constraint violation is less than options.TolCon.
%         options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',tolerFun,'UseParallel', useparallel); 
        options =  optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotstopping},'TolFun',tolerFun,'UseParallel', useparallel); 
        if tProbe.Coil.emfNear == 0
            fitnessfcn = @(x)Final_Tuned_solver(x,tProbe,Visual_glob);
        else
            fitnessfcn = @(x)Final_Tuned_solver_near_inner(x,tProbe,Visual_glob);
        end
        tProbe.Type = 'Tuned';
    else
        options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',tolerFun,'UseParallel', useparallel);
        fitnessfcn = @(x)Final_solver(x,tProbe,Visual_glob);
        tProbe.Type = 'Non-tuned';
    end
    % ,'InitialPopulation',[20,20,0.0004,0.0003]
    % Call the optimization function
    if seeProfile
        profile on
    end
    [x,fval,exitflag,output,population,scores] = ga(fitnessfcn,nvars,A,b,[],[],...
        lb,ub,@simple_nonlinear_constraint,IntCon,options);
    %% Show the results
    if seeProfile
        profile off
        profile viewer
    end
    tProbe.Coil.N_v = x(1); % Loops per layer
    tProbe.Coil.N_l = x(2); % Layers
    if tProbe.Coil.discreteWire == 1 % Increase wire and coating sizes in steps of 0.1 mm
        tProbe.Coil.do = (x(3)+x(4))*1e-4;  % [METERS] wire outer diameter insulation included
        tProbe.Coil.di = x(4)*1e-4;  % [METERS] cooper outer diameter        
    else
        tProbe.Coil.do = x(3)+x(4);  % [METERS] wire outer diameter insulation included
        tProbe.Coil.di = x(4);  % [METERS] cooper outer diameter        
    end

    if tProbe.Coil.emfNear == 0
        tProbe.Coil.Rin = tProbe.Coil.Rout_glob - tProbe.Coil.N_l*tProbe.Coil.do; % Coil internal radius    f0 = (max(freqs_glob)+min(freqs_glob))/2;
    else
        tProbe.Coil.Rin = tProbe.Coil.Rin_glob;
    end
    
    tProbe.Sens.OptProcTime = toc - elapsedTime;
    % #####################################################################################################################################################################
    disp(['Iteration Nr.: ' num2str(countType + ((sweepCount-1)*maxCountType) + ((count-1)*maxSweepCount*maxCountType))])  
    disp('#######################################################################################################################################################')
    if tProbe.Coil.Rin ~= tProbe.Coil.Rin_glob
       warning('3 Main: Rin different to Rin differes from Rin_glob') 
    end
    StrProbeOpt = getProbeDetails(tProbe,graphs);
    if StrProbeOpt.Coil.Rin ~= StrProbeOpt.Coil.Rin_glob
       warning('4 Main: Rin different to Rin differes from Rin_glob') 
    end
    
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
    disp(['Time reference: ' datestr(now)]);
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
if turn_off_PC ==1
    system('shutdown -s -t 180') % To stop it type >> system('shutdown -a')
end

VisualizeCoilOptRes()
%% To visualize one of the results:
% getProbeDetails(ProbeOpt(3,1,1),1); % ProbeOpt(X) where 'X' is the nr. of
% the solution (e.g., 1 for the iteration 1). Adapt acordingly

%% To visualize other parameters by hand:
% [s_test,Probe_test] = CalcSensFromDesign(38,47,0.00029,0.00028,1) ; % (NrLoopsPerLayer,NrLayers,OutWireDiam,CoopDiam, AmplMode)
[s_test,Probe_test] = CalcSensFromDesign(ProbeOpt(3,1,1).Coil.N_v,ProbeOpt(3,1,1).Coil.N_l,ProbeOpt(3,1,1).Coil.do,ProbeOpt(3,1,1).Coil.di,1) ; % (NrLoopsPerLayer,NrLayers,OutWireDiam,CoopDiam, AmplMode)