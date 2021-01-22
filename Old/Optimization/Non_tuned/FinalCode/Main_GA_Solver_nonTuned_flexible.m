% GA optimization 
% 4 inputs 
clear all, close all
fitnessfcn = @Final_solver;
nvars = 4;
% --------INPUTS---------
%     Naxial = x(1);
%     Nradial = x(2);
%     Dtot = x(3);    
%     Dwire = x(4);
%% Linear inequality constraints of the form
% A*x ? b.
A = [0,0,-0.95,1]; %    1*Dwire ? 0.95*Dtot
b = [0];    %0.03 is the value of Rmax
IntCon = [1,2];
% Vector of lower bounds
lb = [1,1,0.0001,0.0001];
% Vector of upper bounds
ub = [50,50,0.01,0.002];
% options =  gaoptimset('PlotFcns',@gaplotbestf);

% delete(gcp)
if exist('poolobj')
    delete(poolobj)
end

poolobj = parpool;

% Define the options
options =  gaoptimset('PlotFcns',@gaplotbestf,'TolFun',1e-16,'UseParallel', true); % The 'TolFun' should be set at least one order 
% of magnitude smaller than the expected target value
% ,'InitialPopulation',[20,20,0.0004,0.0003]
% Call the optimization function
[x,fval,exitflag,output,population,scores] = ga(fitnessfcn,nvars,A,b,[],[],...
    lb,ub,@simple_nonlinear_constraint,IntCon,options);
% [x,fval,exitflag,output,population,scores] = ga(fitnessfcn,nvars,[],[],[],[],...
%     lb,ub,@simple_constraint,IntCon,options)
% kill the parallel computing

delete(poolobj)