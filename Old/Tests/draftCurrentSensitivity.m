% Coil sensitivity acccording to Tashiro 2010, draft

clear all, close all, clc

%% Design by number of turns
di = 0.30e-3; % [METERS] cooper outer diameter
do = 0.33e-3; % [METERS] wire outer diameter insulation included
delta = (do-di)/2; % 0.0001; % [METERS] isolation thickness
a = do/2; % Wire outer radius
N_v = 40; % Number of turns per layer
N_l = 40; % Number of layers
dr = 0.0e-3; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
dz = 0.0e-3; % azimut filling factor
s = do + sqrt(dr^2+dz^2); % distance between wire centres
R_i = 20.07e-3; % Coil internal radius
f = 1/2; %  hexagonal (intercalated)  (f=1)  and orthogonal (f=1/2)


%% Definition of material properties
% Copper
mu_c = 4*pi*1e-7;% relative magnetic permeavility [H/m] of an air core inductor
ro_c = 1.7e-8; % resistivity of the wire

%% Frequencies of interest
freq= logspace(0,9,40);

%% Calculations
[rc,zc] = XY_wire_ort(N_v,N_l,R_i,do,dr,dz);        % Get the location of each wire
figure(1),plot_wires(rc,zc,di,do);                  % Visualize wires for visual confirmation
l_w = coilLength(rc);
R_dc = CoilResistDC(l_w,di);                        % DC Resistance calculation
[R_ac,F,G_2] = CoilResistAC(R_dc,rc,zc,di,s,freq);  % AC Resistance calculation
figure(2),plotResist(R_dc,R_ac,F,G_2,freq);         % Visualize Resistance vs. frequency
[L,Lo,M]= CoilInductanceDC(s,rc,zc);                % Inductance (DC)
[Cs]= CoilCapacitanceDC(R_i,di,do,N_v,N_l,dz,dr);   % Capacitance

% More accurate prediction loop by loop
VonBcoil = 1j*2*(pi^2)*freq*sum(rc.^2);
IonB=(VonBcoil./R_ac).*(1./(1+1j*(2*pi*L*freq./R_ac)));
% % Simplification assuming that all the loops have an average radio 'a'
VonBcoil2 = 1j*2*(pi^2)*freq*N_v*N_l*(mean([max(rc),min(rc)])^2);
IonB2=(VonBcoil2./R_ac).*(1./(1+1j*(2*pi*L*freq./R_ac)));
% Simplified to Rdc
IonB3=(VonBcoil./R_dc).*(1./(1+1j*(2*pi*L*freq./R_dc)));
% Visialize
figure(3), loglog(freq,abs(IonB)), hold on,
loglog(freq,abs(IonB2),'r')
loglog(freq,abs(IonB3),'--g')
grid on, xlabel('Frequency (Hz)'),ylabel('Sensitivity (Ampere/Tesla)')
legend('Accurate sensitivity (loop by loop and R_a_c)','Simplified to an average radius','Simplified to R_d_c')



