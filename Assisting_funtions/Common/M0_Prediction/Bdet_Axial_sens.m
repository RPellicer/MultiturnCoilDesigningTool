%% The fuction for Total magnetization (sensitivity) over Region Of Interest (ROI)

% This is the flux accros the coil (remaining to multiply by the magnetization) when estimating the EMF!!
function [sens_point] = Bdet_Axial_sens(h_coil,A,h_point)
u0 = 4*pi*1e-7;          %magnetic permeability (in S.I. units)
I = 1;                   %current int the coil (Amp)
B_point = u0.*A.^2*I./(2*((h_point+h_coil).^2+A.^2).^(3/2)); % from
sens_point = sum(B_point); % Sum from all the COILS
end