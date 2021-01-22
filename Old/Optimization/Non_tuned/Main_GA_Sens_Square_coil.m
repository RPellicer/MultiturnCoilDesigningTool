%% GA Optimization test 
% for finding
% Number of wires and percentage of copper
% It assumes the shape to be square
clc
clear
close all
fitnessfcn = @Square_solver_far_emf;
nvars = 2;
% use linear inequalities A and b or use Nonlinear by defining function
% Integer constraint
IntCon = [1];
% limits
% 1. Number of wires = 1 to 50
% 2. percent copper = 50% - 90%
lb = [1,40];
ub = [50,90];
% options = psoptimset('TolFun',1e-15);

% Prepare for the parallel computing
if exist('poolobj')
    delete(poolobj)
end
poolobj = parpool;

% Define the options
options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',1e-14,'UseParallel', true);

% Call the optimization function
[x,fval,exitflag,output,population,scores] = ga(fitnessfcn,nvars,[],[],[],[],...
    lb,ub,[],IntCon,options);

% kill the parallel computing
delete(poolobj)