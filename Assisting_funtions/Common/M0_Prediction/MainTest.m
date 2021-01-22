clear all, close all,
I=1;
z_area = 0;
r_area = 0;
R = 0.02;
u0 = 4*pi*1e-7;          %magnetic permeability (in S.I. units)
k2 = abs(4*r_area.*R ./(z_area.^2 +(R+r_area).^2));         %input for elliptic integrals
[K,E] = ellipke(k2);
Bz = ((u0*I)/(2*pi))*(1./sqrt(z_area.^2 +(R+r_area).^2)).*((((R.^2) - (z_area.^2) - (r_area.^2))./(z_area.^2 + (r_area-R).^2)).*E + K);

Bzsimp = u0*I*R^2/(2*(R^2+z_area^2)^(3/2));

r_area = linspace(0,0.019999999,100);
Bz2 = ((u0*I)/(2*pi))*(1./sqrt(z_area.^2 +(R+r_area).^2)).*((((R.^2) - (z_area.^2) - (r_area.^2))./(z_area.^2 + (r_area-R).^2)).*E + K);
plot(r_area,Bz2)
V = 2*pi*trapz(r_area,(Bz2.*r_area)')/(pi*max(r_area)^2); % Average_Bz = doubleCylindricalIntegral(Bz)/Disc_Area;

% % Testing Myers 2007 field strength: Theoreticaly it it 9fT
% Bp = 1;
% Vvoxel = 1e-9;
% z_area = 25e-3;
% R = (65e-3)/2;
% Bdet = u0*2*Vvoxel*netMagnetization(Bp)/(4*pi*(R^2+z_area^2)^(3/2));
% % Bdet = 9.0381e-15, Good!

%% Comparing equations of on axis with off axis 
z_ROI = 0.02;
r_ROI = 0.000001;
R_coil = 0.02;
H_coil = -0.01;
I=1;

Bzon = u0*I*2*pi*(R_coil^2)/(4*pi*(z_ROI^2+R_coil^2)^(3/2));

k2 = abs(4*r_ROI.*R_coil ./(z_ROI.^2 +(R_coil+r_ROI).^2));         %input for elliptic integrals
[K,E] = ellipke(k2);
Bzoff = ((u0*I)/(2*pi))*(1./sqrt(z_ROI.^2 +(R_coil+r_ROI).^2)).*((((R_coil.^2) - (z_ROI.^2) - (r_ROI.^2))./(z_ROI.^2 + (r_ROI-R_coil).^2)).*E + K);

%% Testing the sensitivities to a point or to the volume
z_ROI = 0.02;
r_ROI = 0.000001;
R_coil = 0.02;
H_coil = -0.01;

sens_vol = Bdet_ROI_averag(H_coil,R_coil,r_ROI,z_ROI,z_ROI);
sens_point = Bdet_Axial_sens(H_coil,R_coil,z_ROI);
% sens_vol = Bdet_ROI_averag([H_coil,H_coil],[R_coil,R_coil],r_ROI,z_ROI,z_ROI);
% sens_point = Bdet_Axial_sens([H_coil,H_coil],[R_coil,R_coil],z_ROI);
% They should be equal when the volume is very small

%% Calculating the emf
Bp = 0.1;
z_ROI1 = 10e-3;
z_ROI2 = 30e-3;
r_ROI = 20e-3;
R_coil = 20e-3;
H_coil = 0;
freq = 2e3;

w =2*pi*freq;
sens_vol = Bdet_ROI_averag(H_coil*ones(1,500),R_coil*ones(1,500),r_ROI,z_ROI1,z_ROI2);
sens_point = Bdet_Axial_sens(H_coil*ones(1,500),R_coil*ones(1,500),(z_ROI1+z_ROI2)/2);
% sens_vol = Bdet_ROI_averag(H_coil,R_coil,r_ROI,z_ROI1,z_ROI2);
% sens_point = Bdet_Axial_sens(H_coil,R_coil,(z_ROI1+z_ROI2)/2);
VolumePhant = abs(z_ROI2-z_ROI1)*pi*r_ROI^2;
emf_vol = w*sens_vol*netMagnetization(Bp)*VolumePhant;
emf_point = w*sens_point*netMagnetization(Bp)*VolumePhant;