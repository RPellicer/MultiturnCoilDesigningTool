% Calculating the THERMAL voltage NOISE in series to the inputed resistance in V/sqrt(Hz)
function [en_r] = en_thermal(R)
    kb = 1.38e-23; %k is Boltzmann’s constant (J/K)
    T = 300; % room temperature in kelvin
    en_r = sqrt(4*kb*T*R);
end