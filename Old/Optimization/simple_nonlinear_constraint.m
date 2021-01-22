%% this function is used to apply NON linear constraints
function [c, ceq] = simple_nonlinear_constraint(x)
    setGlobalVariables(); % This function sets the global variables
%    c = [x(2)*x(3)-0.025];       % as 0.03 is the Rmax
    if tProbe.Coil.discreteWire == 1
        c = [x(2)*(x(3)+x(4))*1e-4-tProbe.Coil.Rout_glob + 5e-3;...       % Ensure that the outer diameter of the coil is not bigger than tProbe.Coil.Rout_glob plus a tolerance
            x(1)*(x(3)+x(4))*1e-4-tProbe.Coil.Lmax_glob];       % Ensure the height of the coil is not overpassed.  Ruben 2020 discretise wire choices
    else
        c = [x(2)*x(3)-tProbe.Coil.Rout_glob + 5e-3;...       % Ensure that the outer diameter of the coil is not bigger than  tProbe.Coil.Rout_glob plus a tolerance
            x(1)*x(3)-tProbe.Coil.Lmax_glob];       % Ensure the height of the coil is not overpassed.  
    end
   ceq = [];
end