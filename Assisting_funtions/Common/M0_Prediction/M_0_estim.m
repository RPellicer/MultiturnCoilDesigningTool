%% Magnitude of magnetization of the sample M_0
function [M_0, Bsample] = M_0_estim(Bp,Pd)
    % Field strength Bp in Tesla
% Pd = 6.7e28; of water
    % Other constants
%     h = 6.626e-34; % Js : Planck’s constant
%     I = 1/2; % nuclear spin
%     gyrom = 42.576e6; %  ( Hz T?1) hydrogen's gyromagnetic ratio
%     Ts = 310; % Sample temperature in kelvin
%     k = 1.38e-23; % J/K Boltzmann constant 
%     M_0 = (ones(size(Bp,1),1)*Pd).*((gyrom^2)*(h^2)*I*(I+1)*Bp./(3*k*Ts)); % Amperes/meter
%     mu_0 = 4*pi*1e-7; % = 1.2566e-06 H/m 
% M_0 = (ones(size(Bp,1),1)*Pd).*((42.576e6^2)*(6.626e-34^2)*1/2*(1/2+1)*Bp./(3*1.38e-23*310)); % Amperes/meter
    M_0 = 4.6508e-32*(ones(size(Bp,1),1)*Pd).*Bp; % Amperes/meter
%     Bsample = M_0*mu_0;
    Bsample = M_0*1.2566e-06;    
end
%     n = 6.023e23; % mol^-1 : Avogadro’s number
%     density = 1e6; % g/m^3
%     atomic_mass = 18; % g/mol atomic_mass of water
    % Spin population = Pd water = 6.7e28
%     N = 2*density*n/atomic_mass; % number of spins per unit volume (water has 2 hydrogen molecules and 90% of a human is water)

    % For the protons in water, ? = 6.69 × 10^28 protons/m3. The spin density of human tissue is roughly 0.75 times that of water.
    % Calculation
%     M_0 = N*(gyrom^2)*(h^2)*I*(I+1)*Bp./(3*k*Ts); % Amperes/meter
%   M_0 = B_0*(gy^2)*(h_dash^2)*Pd/(4*k*T); % Equivalent to Mo = B_0*gy^2*h^2*Pd*I*(I+1)/3*k*T for protons (I = 1/2) from ("The NMR receiver")

    %% Ratio M_0/Bp
    % ratio = N*(gyrom^2)*(h^2)*I*(I+1)/(3*k*Ts);
    %  B = mu_o(H+M) [Tesla, T] => the magnetic field created by M is = mu_o*M


%% Tested from another source to double check (Myer et al, JMR 2007)
% ro = 6.7e28; % Spin density
% M_02 = ro*(gyrom^2)*(h^2)*Bp/(4*k*Ts); % Equation from Myer et al, JMR 2007
% Another reference: M0 = (2.34×10-3 A m-1 T-1) Bp. / Busch-2011-thesis-Ultra-low field MRI of prostate cancer using SQUID detection

