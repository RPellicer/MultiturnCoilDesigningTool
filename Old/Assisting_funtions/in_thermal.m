% Calculating the THERMAL current NOISE parallel to the inputed resistance in A/sqrt(Hz)
function [in_r] = in_thermal(R)
    kb = 1.38e-23; %k is Boltzmann’s constant (J/K)
    T = 300; % room temperature in kelvin
    in_r = sqrt(4*kb*T./R);
end