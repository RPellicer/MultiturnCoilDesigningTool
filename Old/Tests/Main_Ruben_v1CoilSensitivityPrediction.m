%% Calculating the sensitivity of the multiloop coil
close all, clear all, clc

%Define global variables
global u0
u0=4*pi*1e-7; %permeability of free space
T = 25 + 273.15; % 25 C in Kelvins
k = 1.3806488e-23; % J/K

%% Calculate the teoretical sensitivity of the sensor
%% Electrical model
% Design by number of turns
di = 0.20e-3; % [METERS] cooper outer diameter
do = 0.24e-3; % [METERS] wire outer diameter insulation included
delta = (do-di)/2; % 0.0001; % [METERS] isolation thickness
a = do/2; % Wire outer radius
N_v = 42; % Number of turns per layer
N_l = 35; % Number of layers
dr = 0.0e-3; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
dz = 0.0e-3; % azimut filling factor
s = do + sqrt(dr^2+dz^2); % distance between wire centres
R_i = 10e-3; % Coil internal radius
f = 1/2; %  hexagonal (intercalated)  (f=1)  and orthogonal (f=1/2)

% Definition of material properties
% Copper
% mu_c = 4*pi*1e-7;% relative magnetic permeavility [H/m] of an air core inductor
% ro_c = 1.7e-8; % resistivity of the wire

% Frequencies of interest
freqtmp2= logspace(0,8,4000);

% Calculations
[rcb,zcb] = XY_wire_ort(N_v,N_l,R_i,do,dr,dz);        % Get the location of each wire
figure(1),plot_wires(rcb,zcb,di,do);                  % Visualize wires for visual confirmation
l_w = coilLength(rcb);
Rb_dc = CoilResistDC(l_w,di);                        % DC Resistance calculation
[Rb,F,G_2] = CoilResistAC(Rb_dc,rcb,zcb,di,s,freqtmp2);  % AC Resistance calculation
figure(2),plotResist(Rb_dc,Rb,F,G_2,freqtmp2);         % Visualize Resistance vs. frequency
[Lb,Lo,M]= CoilInductanceDC(s,rcb,zcb);                % Inductance (DC)
[Cb]= CoilCapacitanceDC(R_i,di,do,N_v,N_l,dz,dr);   % Capacitance

%% Plug here the desired coil parameters unless intented to use the above calculated ones
% For above calculated use
w = 2*pi*freqtmp2;
Rs = Rb;
Ls = Lb;
Cs = Cb;
Seq = sum(pi.*(rcb.^2)); 
% % Or plug here other ones (to check Tashiro's pares for example)
% Rs = 20.*ones(size(freqtmp2));
% Ls = 70e-3;
% Cs = 2.6e-12;
% % My 1st coil
Rs = 71.7.*ones(size(freqtmp2));
Ls = 57.4e-3;
Cs = 8.2e-12;
% % Plain resistors
% Rs = 1000.*ones(size(freqtmp2));
% Ls = 1e-12;
% Cs = 8.2e-12;

% Seq = 1./(1j.*w);

%% Magnetic model
% emf_far=  ?B ?_(i=1:N_t)[?*r_i^2]
emf_far = 1j*w.*Seq; % per B unity!

%% Electric model amplification (Timofeva)
tb = Ls./Rs;
% Amplifier resistors
R1 = 2*5e3; % 2* added for the three OPA Instrumentation amplifier configuiration
R2 = 6e3;
R3 = 6e3;
dV_dB = (emf_far./(Rs.*(1+1j*w.*tb)))*(R1/R2)*R3; % Timofeva 2011, EQ. T_flux without compensation stage
figure, loglog(freqtmp2,abs(dV_dB),'b'); 
grid on, title('Field amplification Vs. frequency')

%% Add noise model here
% Amplifier specs. The noise model is taken from "Calculating noise figure in op amps" (Texas Instrument) which matches with 
% *The circuits and filters handbook* FIGURE 15.27 in Chapter 15
i_n = 0.8e-12; %INA217
e_n = 1.2e-9;  %INA217
e_OutStagNoise = 90e-9;
% i_n = 1e-12;
% e_n = 0.9e-9;
Zsen = (Rs + 1j*w.*Ls)./((1-Ls*Cs*w.^2)+1j*w.*Rs*Cs); % Z_Thevenin Acording to Timofeeva et al, 2011 -"Differential Search Coils Based Magnetometers: Conditioning, Magnetic Sensitivity, Spatial Resolution"
Zs = Rs + 1j*w.*Ls; % Impedance of the branch (inductance + resistance). The other branch (i.e., the inductance) doesn't account in this
Ys = 1./Zsen;
Zf = R1;
Rf = R1;
Yf = 1./Zf;
% Define amplifier gains
% sigGain = abs(Zf./Zs);
sigGain = abs(Zf./Zs);
noiseGain = abs(1+(Zf./Zsen));
% % Calculate individually the VOLTAGE noise sources INPUT-REFERRED
% ce_ni = abs(e_n.*(Zs./Zf+1));                  % Voltage noise of the amplifier
% ci_ii = abs(i_n.*Zs);                         % Current noise of the amplifier
% ce_s = abs(sqrt(4*k*T.*Rs));            % Thermal noise of the source
% ce_f = abs(sqrt(4*k*T*Rf).*Zs./Zf);      % Thermal noise of the feedback impedance
% en_out = ce_ni .* sigGain;
% ei_out = ci_ii .* sigGain;
% es_out = ce_s .* sigGain;
% ef_out = ce_f .* sigGain;

% ce_ni = e_n.*((Zs./Zf)+1);                  % Voltage noise of the amplifier
ce_ni = e_n.*noiseGain./sigGain;                  % Voltage noise of the amplifier
ci_ii = i_n.*Zs;                         % Current noise of the amplifier
% ce_s = sqrt(4*k*T.*Rs);            % Thermal noise of the source
ce_s = en_thermal(Rs);
ce_f = en_thermal(Rf+(4*6e3)).*abs(Zs./Zf);     % Thermal noise of the    
                                                % feedback impedance. In
                                                % the case if the INA217 it
                                                % is 2x5k ohm + 4x6k ohm
% en_out = abs(ce_ni .* sigGain);
en_out = e_n.*noiseGain;
ei_out = abs(ci_ii .* sigGain);
es_out = abs(ce_s .* sigGain);
ef_out = abs(ce_f .* sigGain);
oe = e_OutStagNoise.*ones(size(freqtmp2)) ;  % This noise dominates with small gains                                              

% ce_in = sqrt(ce_ni^2 + ci_ii.^2 + ce_s.^2 + ce_f.^2); % VOLTAGE noise sources INPUT-REFERRED
vn_out = sqrt(en_out.^2 + ei_out.^2 + es_out.^2 + ef_out.^2 + oe.^2);   % The noise is calculated to the output with the 
                            % signal gain because the noise gain was taken into account for the calculation 
                            % of each of the sources as in "Calculating noise figure in op amps".

% Comparing with Netzer et al. 1982, "The Design of Low-Noise Amplifiers"
currGain = abs(Zf).*ones(size(freqtmp2));
% i_in2 = sqrt(i_n^2 + e_n^2.*(Ys+Yf).^2 + 4*k*T.*real(Ys+Yf)); % this is the original equation
en_out2 = abs(e_n.*abs(Ys+Yf)).*currGain;
ei_out2 = abs(i_n).*currGain;
esf_out2 = abs(sqrt(4*k*T.*real(Ys+Yf))).*currGain;
es_out2 = abs(sqrt(4*k*T.*real(Ys))).*currGain;
% es_out2 = abs(sqrt(4*k*T.*1./Rs)).*currGain;
ef_out2 = abs(sqrt(4*k*T.*real(Yf))).*currGain;
vn_out2 = sqrt(en_out2.^2 + ei_out2.^2 + es_out2.^2 + ef_out2.^2 + ef_out2.^2);
% vn_out2 = i_in2 * sigGain;

%Comparing it with Tashiro (2006)"Optimal design of an air-core..." Eqs 10-14
noiseGain3 = (1 + Rf./Rs);
En3 = e_n;
Ei3 = abs(i_n*((Rs.*Rf)./(Rs+Rf)));
Er3 = abs(sqrt(4*k*T*((Rs.*Rf)./(Rs+Rf))));
% vn_out3 = noiseGain*sqrt(En3.^2 + Ei3.^2 + Er3.^2);
en_out3 = En3.*noiseGain3;
ei_out3 = Ei3.*noiseGain3;
es_out3 = Er3.*noiseGain3;
vn_out3 = sqrt(en_out3.^2 + ei_out3.^2 + es_out3.^2);
% Visualize gains
figure, loglog(freqtmp2,abs(sigGain),'b');
hold on, grid on,
loglog(freqtmp2,noiseGain,'r');
loglog(freqtmp2,noiseGain3,'g');
legend('Signal gain','Noise gain','Noise gain Tashiro');
title('Compare gains')
% Compare noise models
figure, loglog(freqtmp2,abs(vn_out),'b');
hold on, grid on
loglog(freqtmp2,abs(vn_out2),'-.r');
loglog(freqtmp2,abs(vn_out3),'--k');
legend('Texas Instrument noise model','Netzer noise model','Tashiro noise model')

% Compare each component
figure, 
loglog(freqtmp2,abs(vn_out),'-*k');
hold on, grid on
loglog(freqtmp2,abs(en_out),'m');
loglog(freqtmp2,abs(ei_out),'r');
loglog(freqtmp2,abs(es_out),'b');
loglog(freqtmp2,abs(ef_out),'g');
loglog(freqtmp2,abs(oe),'c');
legend('Total output noise','en_out','ei_out','es_out','ef_out','OutputStageNoise')
title('Texas Instrum.')
figure, 
loglog(freqtmp2,abs(vn_out2),'-*k');
hold on, grid on
loglog(freqtmp2,abs(en_out2),'m');
loglog(freqtmp2,abs(ei_out2),'r');
loglog(freqtmp2,abs(es_out2),'b');
loglog(freqtmp2,abs(ef_out2),'g');
legend('Total output noise','en_out','ei_out','es_out','ef_out'),
title('Netzer')
figure, 
loglog(freqtmp2,abs(vn_out3),'-*k');
hold on, grid on
loglog(freqtmp2,abs(en_out3),'m');
loglog(freqtmp2,abs(ei_out3),'r');
loglog(freqtmp2,abs(es_out3),'b');
legend('Total output noise','en_out','ei_out','es_out'),
title('Tashiro')

% Compare gains
figure, loglog(freqtmp2,abs(Zs))
hold on, grid on,
loglog(freqtmp2,real(Zs))
loglog(freqtmp2,Rs)
loglog(freqtmp2,abs(sigGain))
loglog(freqtmp2,noiseGain)
title('Gains and impedance visualization')
legend('Zs','Real(Zs)','Rs','Signal gain','Noise gain')

% Compare field gain with noise floor
figure, loglog(freqtmp2,abs(dV_dB),'b'); 
hold on, grid on, title('Field amplification Vs. frequency and noise floor')
loglog(freqtmp2,abs(emf_far.*sigGain))
loglog(freqtmp2,vn_out)
legend('Field amplification Timofeeva','Field amplification Ruben','Noise floor Texas Instruments')
figure,loglog(freqtmp2,abs(vn_out./(emf_far.*sigGain)))
grid on, hold on, title('Sensitivity Vs. Frequency')
ylabel('Sensitivity (T/Hz^1^/^2)'), xlabel('Frequency (Hz)')

% realistic sensitivity Figure for Michael
load('temp1')
% G2temp = interp1(freqtmp2,abs(emf_far.*sigGain),f2);
% loglog(dVdTtina(:,1),dVdTtina(:,2));
G2temp = interp1(dVdTtina(:,1)',dVdTtina(:,2)',f2);
h=figure;
loglog(f2,G2temp,'r')
hold on
loglog(freqtmp2,abs(emf_far.*sigGain*100),'--g');% The "*100" comes from the gain of the second stage
title('comparison of gains with and without filters(from TINA software)'); %  TNA file "Ruben_INA217_I-V" that includes the filters with the gain without filters

% Sens = abs(sqrt(Pxx2)./(G2temp*100)); 
Sens = abs(sqrt(Pxx2)./(G2temp)); 
figure,
loglog(f2,Sens,'r');
hold on
loglog(freqtmp2,abs(vn_out./(emf_far.*sigGain)),'b')
grid on
legend('Measured','Mean value')
temp4 = f2>=1500 & f2<=30000; % This sets the range we want to visualize chopping off the effects of the filters
indtemp = find(temp4);
f3 = f2(indtemp);
G3temp = G2temp(indtemp);
measSig3 = abs(sqrt(Pxx2(indtemp)));
loglog(f3,(measSig3./G3temp),'g');

% Arrange a nice figure for the grant
figure,
loglog(freqtmp2,abs(vn_out./(emf_far.*sigGain)),'b')
hold on
grid on
loglog(f3,(measSig3./G3temp),'g');
legend('Simulated','Measured'),title('Sensitivity Vs. Frequency of the multiturn sensor')
ylabel('Sensitivity (T/Hz^1^/^2)'), xlabel('Frequency (Hz)')

