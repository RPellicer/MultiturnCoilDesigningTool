%% This script helps calculating multiple turn surface coil characteristics
% From Martinez et all, 2014 “On evaluation of inductance, dc resistance and capacitance of coaxial inductors at low frequencies”
clear all, close all, clc

%% Design by number of turns
d_0 = 0.25e-3; % [METERS] cooper outer diameter
delta = 0.0e-3;% 0.0001; % [METERS] isolation thickness
a = d_0/2 + delta; % Wire outer radius
N_v = 20; % Number of turns per layer
N_l = 20; % Number of layers
ksi_r = 0.0e-3; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
ksi_z = 0.0e-3; % azimut filling factor
R_i = 10e-3; % Coil internal radius
f = 1/2; %  hexagonal (intercalated)  (f=1)  and orthogonal (f=1/2)

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
tic
for i = 1:N_l
    for j = 1:N_v
            pos(i,j,1) = R_i+(2*i-1)*a-(i-1)*ksi_r;
            pos(i,j,2) = (j-1/2)*(2*a+ksi_z);
    end
end
if f==1;
%     pos(:,:,1) = (pos(:,:,1)-(R_i+a))*sin(pi/3)+(R_i+a); % Brings the wire to be in contact
    for i = 2:2:N_l
        pos(i,:,2) = pos(i,:,2)+a+ksi_z/2;
    end
end
toc
posvect = [reshape( pos(:,:,1).' ,1,numel(pos(:,:,1))) ; reshape( pos(:,:,2).' ,1,numel(pos(:,:,2)))]';
circles(squeeze(ones(size(posvect,1),1)*(d_0/2)), [posvect(:,2),posvect(:,1)],'color','b');
hold on
plot(0,0,'*r'), xlabel(['Z (m) N_v=' num2str(N_v)]), ylabel(['r (m) N_l=' num2str(N_l)])
axis equal, grid on
hold off
toc
%% DC Resistance calculation
% the next simplified equation is not too accurate so we use the more
% complete version
% l_w = pi*N_t*(2*R_i+(2*N_l-1)*a-(N_l-1)*ksi_r-(N_v-1)*ksi_z); % Wire length, eq. 16, Martinez et all 2014
toc
l_w = 0;
for i = 1: N_l
    for j = 1: N_v
        l_w = l_w + 2*pi*(R_i + (2*i - 1)*a-(i-1)*ksi_r-(j-1)*ksi_z);
    end
end
R_dc = (4*ro_c/(d_0^2))*(l_w/pi); % Total DC resistance of the cable
toc
%% Total inductance of the coil
% Position  of  each  turn  within  the  array
% for i = 1:N_t
%     r(i) = R_i + (2*i-1)*a - (i-1)*ksi_r;
% end
% r = repmat(r,[1 N_v])'; % Other turns will be as far as this first calculation
k = 4*(posvect(:,1).*(posvect(:,1)-a))./((2.*posvect(:,1)-a).^2);
% k = 4*(r.*(r-a)./((2*r-a).^2));
%Complete elliptic integrals of first and second kind
[K,E] = ellipke(k);
% Self-inductance (Independent inductance of each turn)
L_0mem = [];
L_0memWiki = [];
Y = 0; % Y is a constant. Y=0 when the current flows in the surface of the wire (skin effect), Y=1/4 when the current is homogeneous across the wire.
for i=1:N_t
    L_0mem(i) = mu_c*(2*posvect(i,1)-a)*((1-(1/2)*(k(i)^2))*K(i)-E(i));
    L_0memWiki(i) = mu_c.*(posvect(i,1)).*(log(8.*posvect(i,1)./a)-2);
end
L_0 = sum(sum(L_0mem));
L_0Wiki = sum(sum(L_0memWiki));
% Mutual inductance

km = [];
for i = 1:N_t
    for j = 1:N_t
        if i~=j
            km(i,j) = (4*posvect(i,1)*posvect(j,1))/((posvect(i,1)+posvect(j,1))^2+(posvect(i,2)-posvect(j,2))^2);
        end
    end
end

[Km,Em] = ellipke(km);

% Mmem = [];
% M = 0;
% for i = 1:N_t
%     for j = 1:N_t
%         if i~=j 
%             Mmem(i,j) = 2*mu_c*(sqrt(posvect(i,1)*posvect(j,1))/km(i,j))*((1-(1/2)*(km(i,j)^2))*Km(i,j)-Em(i,j));
%         end
%     end
% end
% M = sum(sum(Mmem));


% % ################## In matrix ways as Martinez et all #################
% 
Km_NTxNT = (1-((km.^2)./2)).*Km;
Em_NTxNT = Em;
Fm_NTxNT = Km_NTxNT-Em_NTxNT;
Ones_NTx1 = ones(N_t,1);
% Mm_NTxNT = [];
% for i = 1:N_t
%     for j = 1:N_t
%         if i~=j 
%             Mm_NTxNT(i,j) = sqrt(posvect(i,1)*posvect(j,1))/km(i,j);
%         end
%     end
% end
Mm_NTxNT = sqrt(posvect(:,1)*posvect(:,1)')./km;
for i = 1 :N_t
    Mm_NTxNT(i,i) = 0;
end

% M_Martinez = 2*mu_c*((Mm_NTxNT.*Fm_NTxNT)*Ones_NTx1)'*Ones_NTx1;
M_Martinez = 2*mu_c*(Mm_NTxNT.*Fm_NTxNT);
M = (M_Martinez*Ones_NTx1)'*Ones_NTx1;
toc
% 
% % ######################################################################
% imshow(Mmem,[0 max(max(Mmem))])
L = L_0 + M;
text = ['L_0 = ' num2str(L_0) '  and  M = ' num2str(M), ' , thus L = ' num2str(L)];
disp(text);