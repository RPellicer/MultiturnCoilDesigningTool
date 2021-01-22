%% Calculation of coil AC characteristics and sensitivity
% from I.M. Savukov et al. / Journal of Magnetic Resonance 185 (2007)
% Detection of NMR signals with a radio-frequency atomic magnetometer
clear all, close all, clc

%% Definition of frequency
freq= logspace(0.1,9,40);
omega = 2*pi*freq;
%% Design by number of turns
d_0 = 0.511e-3; % [METERS] cooper outer diameter
delta = 0.164e-3;% 0.0001; % [METERS] isolation thickness
a = d_0/2 + delta; % Wire outer radius
N_v = 24; % Number of turns per layer
N_l = 60; % Number of layers
ksi_r = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
ksi_z = 0.0; % azimut filling factor
R_i = (20e-3)/2; % Coil internal radius
f = 1/2; %  hexagonal  (f=1)  and orthogonal (f=1/2)
W_1 = (a+ksi_r) * N_l;
W_2 = (a+ksi_z) * N_v;
D = R_i + W_1/2;

%% Definition of material properties
% Copper
mu_c = 4*pi*1e-7;% relative magnetic permeavility [H/m] of an air core inductor
epsilon_r = 1; % relative electrical permittivity;
ro_c = 1.7e-8; % resistivity of the wire
% Isolator
epsilon_ef = 1; % dielectric constant of foil between layers;

% calculation of design dimentions
N_t = N_v * N_l;% total number of turns
%% Position of each wire
% for k = 1:N_l
%     for j = 1:N_v
%             pos(k,j,1) = R_i+(2*k-1)*a-(k-1)*ksi_r;
%             pos(k,j,2) = (j-1/2)*(2*a+ksi_z);
%     end
% end
temp_l = [1:N_l]'*ones(1,N_v);
temp_v = ones(N_l,1)*[1:N_v];
pos(:,:,1) = R_i+(2*temp_l-1)*a-(temp_l-1)*ksi_r; 
pos(:,:,2) = (temp_v-1/2)*(2*a+ksi_z);

if f==1;
    pos(:,:,1) = (pos(:,:,1)-(R_i+a))*sin(pi/3)+(R_i+a);
%     for k = 2:2:N_l
%         pos(k,:,2) = pos(k,:,2)+a+ksi_z/2;
%     end
    k = 2:2:N_l;
    pos(k,:,2) = pos(k,:,2)+a+ksi_z/2;
end

posvect = [reshape( pos(:,:,1).' ,1,numel(pos(:,:,1))) ; reshape( pos(:,:,2).' ,1,numel(pos(:,:,2)))]';
circles(squeeze(ones(size(posvect,1),1)*(d_0/2)), [posvect(:,2),posvect(:,1)],'color','b');
hold on
circles(squeeze(ones(size(posvect,1),1)*(a)), [posvect(:,2),posvect(:,1)],'color','g');
plot(posvect(:,2),posvect(:,1),'+r');
hold on
plot(0,0,'*r'), xlabel(['Z (m) N_v=' num2str(N_v)]), ylabel(['r (m) N_l=' num2str(N_l)])
axis equal, grid on
hold off

%% DC Resistance calculation
% the next simplified equation is not too accurate so we use the more
% complete version
% l_w = pi*N_t*(2*R_i+(2*N_l-1)*a-(N_l-1)*ksi_r-(N_v-1)*ksi_z); % Wire length, eq. 16, Martinez et all 2014
% l_w = 0;
% for k = 1: N_l
%     for j = 1: N_v
%         l_w = l_w + 2*pi*(R_i + (2*k - 1)*a-(k-1)*ksi_r-(j-1)*ksi_z);
%     end
% end
temp_lw = 2*pi*(R_i + (2*temp_l - 1)*a-(temp_l-1)*ksi_r-(temp_v-1)*ksi_z); 
l_w = sum(sum(temp_lw));
R_dc = ro_c*l_w/(pi*(d_0/2)^2); % Total DC resistance of the cable


%% Calculation of the AC resistance
s = 2*a; % s is the spacing between the centers of the wires
posvect_s = posvect/s; % because u is calculated in units of s.
posvect_s(:,1) = posvect_s(:,1)-posvect_s(1,1)+1;
posvect_s(:,2) = posvect_s(:,2)-posvect_s(1,2)+1;
posvect_s = round(posvect_s);

% skin depth
skinDepth = sqrt(2*ro_c./(omega.*mu_c));
z = d_0./(skinDepth*sqrt(2));

% Skin depth EFFECT
F = (-z.^2./8).*imag(besselj(3,z*sqrt(-1i))./besselj(1,z*sqrt(-1i)));

% Proximity EFFECT
G = (-z.^2./8).*imag(besselj(2,z*sqrt(-1i))./besselj(0,z*sqrt(-1i)));
% u = 0;
% for k = 1:N_t
%     temp_u1 = 0;
%     temp_u2 = 0;
%     for j = 1:N_t
%         if k~=j 
%             temp_u1 = temp_u1 + (posvect_s(j,1)-posvect_s(k,1))/((posvect_s(j,1)-posvect_s(k,1))^2+(posvect_s(j,2)-posvect_s(k,2))^2);
%             temp_u2 = temp_u2 + (posvect_s(j,2)-posvect_s(k,2))/((posvect_s(j,1)-posvect_s(k,1))^2+(posvect_s(j,2)-posvect_s(k,2))^2);
%             tempi3(j,k) = (posvect_s(j,1)-posvect_s(k,1))/((posvect_s(j,1)-posvect_s(k,1))^2+(posvect_s(j,2)-posvect_s(k,2))^2);
%             tempi4(j,k) = (posvect_s(j,2)-posvect_s(k,2))/((posvect_s(j,1)-posvect_s(k,1))^2+(posvect_s(j,2)-posvect_s(k,2))^2);
%         end
%     end
%     u = u + (1/N_t)*(temp_u1^2+temp_u2^2);
% end
% 
temp_ul = ones(N_t,1)*posvect_s(:,1)';
temp_uv = ones(N_t,1)*posvect_s(:,2)';
mask = ~diag(ones(N_t,1));
Au = (temp_ul - temp_ul')./((temp_ul - temp_ul').^2+(temp_uv - temp_uv').^2);
Bu = (temp_uv - temp_uv')./((temp_ul - temp_ul').^2+(temp_uv - temp_uv').^2);
u = (1/N_t)*(nansum(nansum(mask.*Au).^2)+...
    nansum(nansum(mask.*Bu).^2));

G_2 = u*((d_0^2)/(s^2))*G;

% Total AC resistance
R_ac = R_dc*(1+F+G_2);
figure(2),
subplot(2,1,1)
semilogx(freq,R_ac,'-*'), ylabel('resistance(ohm)'),xlabel('frequency(Hz)')
hold on
semilogx(freq,F*R_dc,'r')
semilogx(freq,G_2*R_dc,'g')
legend('R_a_c','Skin depth EFFECT','Proximity EFFECT')
hold off
subplot(2,1,2)
loglog(freq,R_ac,'-*'), ylabel('resistance(ohm)'),xlabel('frequency(Hz)')
hold on
loglog(freq,F*R_dc,'r')
loglog(freq,G_2*R_dc,'g')
legend('R_a_c','Skin depth EFFECT','Proximity EFFECT')
grid on
hold off

% Calculated following Sakurov 2007 formula
% Johnson noise of coil R_ac
T = 25 + 273.15; % 25 C in Kelvins
kb = 1.3806488e-23; % J/K
v_Rac=sqrt(4*kb*T*R_ac); % V/Hz^1/2
figure(3)
loglog(freq,v_Rac,'b')
title('Johnson noise of coil'), ylabel('resistive noise(v/Hz^2)'),xlabel('frequency(Hz)')
grid on

% Sensitivity
ro_ACfill = R_ac*W_1*W_2/(l_w*N_t);
V_w = pi*D*W_1*W_2; % coil 3D volume
coilSens = (8./(omega*D)).*sqrt(kb*T*ro_ACfill./V_w);% calculated 4.5fT/sqrt(Hz)~20fT/sqrt(Hz)reported
figure(4)
loglog(freq,coilSens,'b')
title('Sensitivity'), ylabel('Sensitivity(T/Hz^2)'),xlabel('frequency(Hz)')
grid on

% D = 2*(R_i + (a*N_l)); % diamenter in between coil crossection centers
% figure(1), hold on
% plot(0,D/2,'+k')
% hold off
% v_Rac2=sqrt(16*kb*T*ro_c*D*N_t / (d_0^2)); % It is in good agreement 

% % Voltage induced in the coil
% B = 1;
% v = B*omega*pi*(D^2)*N_t/4;

% Sensitivity
% ro_w = v_Rac*
% sensitivity = (8/(omega*D))*sqrt(kb*T*ro/(pi*D))

% test = decimate(R_ac,2);
down_sampling_factor = 4;
R_ac_simpl = round(downsample(R_ac,down_sampling_factor));
freq_simpl = round(downsample(freq,down_sampling_factor));
loglog(freq,R_ac,'*')
hold on
loglog(freq_simpl,R_ac_simpl,'+r')
grid on
phase_simpl = zeros(size(freq_simpl));
interf = [freq_simpl;1./R_ac_simpl;phase_simpl];
interf2 = reshape(interf,1,size(interf,1)*size(interf,2));
table = [];
for i = 1:length(freq_simpl)
    table = [table '(' num2str(interf(1,i)) ',' num2str(interf(2,i)) ',' num2str(interf(3,i)) ')'];
end
    