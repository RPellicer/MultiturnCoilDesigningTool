%% From Calculated signal-to-noise ratio of MRI detected with SQUIDs
% and Faraday detectors in fields from 10 lT to 1.5 T

function M_0 = netMagnetization(B_0)
    h_dash = 6.62606957e-34/(2*pi); % Plank's constant (m^2*kg/s), h_dash=h/(2*pi)!!
    k = 1.38064852e-23; % Boltzmann's constant (m^2*kg/(s^2*K))
    T = 310;            % Temperature in Kelvin (0C = 273.15K) (36.85C = 310K) (24.85C=298K)
    Pd = 6.7e28;        % Proton nuclear density of water (m^3) (The spin density of human tissue is roughly 0.75 times that of water, from  "Ultra-low field MRI of prostate cancer using SQUID detection-thesis-pg68")
    gy = 42.58e6*2*pi;       % Gyromagnetic ratio (T^-1 s^-1)=> gamma_bar = gamma*2*pi = 42.5774806(10) MHz·T?1!!!
    M_0 = B_0*(gy^2)*(h_dash^2)*Pd/(4*k*T); % Equivalent to Mo = B_0*gy^2*h^2*Pd*I*(I+1)/3*k*T for protons (I = 1/2) from ("The NMR receiver")
end