%% Calculation of coil AC characteristics and sensitivity
% from I.M. Savukov et al. / Journal of Magnetic Resonance 185 (2007)
% Detection of NMR signals with a radio-frequency atomic magnetometer
clear all, close all, clc

%% Definition of frequency
freq= logspace(3,9,100);
omega = 2*pi*freq;
%% Design by space as boundary condition
D_in = 10e-3;
D_out = 30e-3;
D = (D_out+D_in)/2;
W_1 = (D_out-D_in)/2;
W_2 = 12e-3;

ksi_r = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
ksi_z = 0.0; % azimut filling factor
R_i = D_in/2; % Coil internal radius
f = 1/2; %  hexagonal  (f=1)  and orthogonal (f=1/2)

% Amplifier
NoiseAmp = 1.2e-9; % 1.2nV/sqrt(Hz) is the noise floor of the pre-amplifier
% Coil properties
% diameter and insulation of interest will be compared
% ####################################################
guessWire = 0.2e-3; % meters
guessInsulation =  0.03e-3; % meters
% ####################################################
% Array of different diameter and insulation posibles
d_0M = logspace(-4,-3,20);% [METERS] cooper outer diameter
deltaM =logspace(-5,-4,20); %  [METERS] insulation thickness
% find the closest values in the arrays to the ones of interest
[~,idx]= min(abs(d_0M-guessWire));
choiceWire  = idx;
[~,idx]= min(abs(deltaM-guessInsulation));
choiceInsulation = idx;
% from array to matrix to have a convination of both different diameter and
% Insulations
d_0M = ones(size(deltaM,2),1)*d_0M;
deltaM = deltaM'*ones(1,size(d_0M,2));

aM = d_0M/2 + deltaM; % Wire outer radius

N_vM = ceil(W_2./(2*aM)); % Number of turns per layer
N_lM = ceil(W_1./(2*aM)); % Number of layers
N_tM = N_lM .* N_vM; % Total number of turns

%% Definition of material properties
% Copper
mu_c = 4*pi*1e-7;% relative magnetic permeavility [H/m] of an air core inductor
epsilon_r = 1; % relative electrical permittivity;
ro_c = 1.7e-8; % resistivity of the wire
% Isolator
epsilon_ef = 1; % dielectric constant of foil between layers;

% Self-resonance turns limit from Savukov 2007 (very much simplified)
c = 299792458; % m/s, Speed of light
turnMaxResFreq = c./(D*omega);

%% Loop calculating the characterization for all the convinations.
if ~matlabpool('size')
    poolobj = parpool(3);
end
bar = ProgressBar(size(deltaM,2));
for o=1:size(deltaM,1)
    parfor p=1:size(d_0M,2)
        d_0 = d_0M(o,p);
        a = aM(o,p);
        N_v = N_vM(o,p);
        N_l = N_lM(o,p);
        N_t = N_tM(o,p);
        %% Position of each wire
        temp_l = [1:N_l]'*ones(1,N_v);
        temp_v = ones(N_l,1)*[1:N_v];
        pos = [];
        pos(:,:,1) = R_i+(2*temp_l-1)*a-(temp_l-1)*ksi_r; 
        pos(:,:,2) = (temp_v-1/2)*(2*a+ksi_z);

        if f==1;
            pos(:,:,1) = (pos(:,:,1)-(R_i+a))*sin(pi/3)+(R_i+a);
            k = 2:2:N_l;
            pos(k,:,2) = pos(k,:,2)+a+ksi_z/2;
        end

        posvect = [reshape( pos(:,:,1).' ,1,numel(pos(:,:,1))) ; reshape( pos(:,:,2).' ,1,numel(pos(:,:,2)))]';

        %% DC Resistance calculation
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

        % Calculated following Sakurov 2007 formula
        % Johnson noise of coil R_ac
        T = 25 + 273.15; % 25 C in Kelvins
        kb = 1.3806488e-23; % J/K
        v_Rac=sqrt(4*kb*T*R_ac); % V/Hz^1/2
        
        % Voltage induced in the coil per testla
        v_ind = omega*pi*D^2*N_t/4;% following Sakurov 2007 formula adapted to V/B

        % Sensitivity
        ro_ACfill = R_ac*W_1*W_2/(l_w*N_t);
        V_w = pi*D*W_1*W_2; % coil 3D volume
        coilSens = (8./(omega*D)).*sqrt(kb*T*ro_ACfill./V_w);% calculated 4.5fT/sqrt(Hz)~20fT/sqrt(Hz)reported

        % Saving the values
        coilSensM(o,p,:) = coilSens; 
        v_RacM(o,p,:) = v_Rac;
        R_acM(o,p,:) = R_ac;
        v_indM(o,p,:) = v_ind;
        zM(o,p,:) = z;
    end
    bar.progress; % Update progress bar
end
bar.stop; % Stop progress bar
save('AC_Coil_Optimization_from_winding_Volume_workspace');

% Plot the sensitivity
h=figure;
test = squeeze(coilSensM(:,choiceWire,:));
loglog(freq,test), grid on
mesleg = [repmat(['Nt = '],size(N_tM,1),1) num2str(N_tM(:,choiceWire))...
    repmat([', isol. thick.= '],size(N_tM,1),1) num2str(deltaM(:,choiceWire)*1000) repmat([' mm'],size(N_tM,1),1)];
legend(mesleg)
title(['Coil sensitivity for different insulation thickness (wire diam.= ' num2str(d_0M(1,choiceWire)*1000) ' mm)'])
xlabel('Freq. Hz'), ylabel('Coil sensitivity T/Hz^1^/^2')
h=figure;
test = squeeze(coilSensM(choiceInsulation,:,:));
loglog(freq,test'), grid on
mesleg = [repmat(['Nt = '],size(N_tM,2),1) num2str(N_tM(choiceInsulation,:)')...
    repmat([', wire diam = '],size(N_tM,2),1) num2str(d_0M(choiceInsulation,:)'*1000) repmat([' mm'],size(N_tM,2),1)];
legend(mesleg)
title(['Coil sensitivity for different wire thickness (isol.= ' num2str(deltaM(choiceInsulation,1)*1000) ' mm)'])
xlabel('Freq. Hz'), ylabel('Coil sensitivity T/Hz^1^/^2')

% Plot the resistance
h=figure;
test = squeeze(R_acM(:,choiceWire,:));
loglog(freq,test), grid on
mesleg = [repmat(['Nt = '],size(N_tM,1),1) num2str(N_tM(:,choiceWire))...
    repmat([', isol. thick.= '],size(N_tM,1),1) num2str(deltaM(:,choiceWire)*1000) repmat([' mm'],size(N_tM,1),1)];
legend(mesleg)
title(['AC resistance  for different insulation thickness (wire diam.= ' num2str(d_0M(1,choiceWire)*1000) ' mm)'])
xlabel('Freq. Hz'), ylabel('Ohms')
h=figure;
test = squeeze(R_acM(choiceInsulation,:,:));
loglog(freq,test'), grid on
mesleg = [repmat(['Nt = '],size(N_tM,2),1) num2str(N_tM(choiceInsulation,:)')...
    repmat([', wire diam = '],size(N_tM,2),1) num2str(d_0M(choiceInsulation,:)'*1000) repmat([' mm'],size(N_tM,2),1)];
legend(mesleg)
title(['AC resistance for different wire thickness (isol.= ' num2str(deltaM(choiceInsulation,1)*1000) ' mm)'])
xlabel('Freq. Hz'), ylabel('Ohms')

% Plot the thermal coil noise (v/sqrt(Hz))
h=figure;
test = squeeze(v_RacM(:,choiceWire,:));
loglog(freq,test), grid on
hold on
line([freq(1) freq(end)],[NoiseAmp NoiseAmp],'Color','r','LineStyle','--','LineWidth',2)
hold off
mesleg = [repmat(['Nt = '],size(N_tM,1),1) num2str(N_tM(:,choiceWire))...
    repmat([', isol. thick.= '],size(N_tM,1),1) num2str(deltaM(:,choiceWire)*1000) repmat([' mm'],size(N_tM,1),1)];
mesleg = [cellstr(mesleg); 'Pre-amplifier`s Noise floor'];
legend(mesleg)
title(['Coil thermal noise for different insulation thickness (wire diam.= ' num2str(d_0M(1,choiceWire)*1000) ' mm)'])
xlabel('Freq. Hz'), ylabel('Coil thermal noise V/Hz^1^/^2')

h=figure;
test = squeeze(v_RacM(choiceInsulation,:,:));
loglog(freq,test'), grid on
hold on
line([freq(1) freq(end)],[NoiseAmp NoiseAmp],'Color','r','LineStyle','--','LineWidth',2)
hold off
mesleg = [repmat(['Nt = '],size(N_tM,2),1) num2str(N_tM(choiceInsulation,:)')...
    repmat([', wire diam = '],size(N_tM,2),1) num2str(d_0M(choiceInsulation,:)'*1000) repmat([' mm'],size(N_tM,2),1)];
mesleg = [cellstr(mesleg); 'Pre-amplifier`s Noise floor'];
legend(mesleg);
title(['Coil thermal noise for different wire thickness (isol.= ' num2str(deltaM(choiceInsulation,1)*1000) ' mm)'])
xlabel('Freq. Hz'), ylabel('Coil thermal noise V/Hz^1^/^2')

matlabpool close

%#########################################################################################################
%% Search for optimal values
w_temp = v_indM./v_RacM; 
%% Detect the setups where amplifier noise is dominant
domin_tresshold = v_RacM<NoiseAmp;
domin_v_indM = v_indM.*domin_tresshold;
w_temp2 = domin_v_indM./v_RacM; 
% The freq nr. 6 is about 2kHz
weight = w_temp(:,:,6);
weight2 = w_temp2(:,:,6); 
h=figure;
hSurface1 = surf(N_tM,weight,deltaM);
set(hSurface1,'FaceColor',[1 0 0],'FaceAlpha',0.5);
view(2)
set(gca,'xscale','log')
set(gca,'yscale','log')
hold on
hSurface = surf(N_tM,weight2,deltaM);
set(hSurface,'FaceColor',[1 0.7 0.7],'FaceAlpha',0.5);
hLine1_10 = scatter3(N_tM(1,1),weight(1,1),deltaM(1,1),'filled','sk');
hLine1_11 = scatter3(N_tM(end,1),weight(end,1),deltaM(end,1),'*k');
scatter3(N_tM(choiceInsulation,choiceWire),weight(choiceInsulation,choiceWire),deltaM(choiceInsulation,choiceWire),'filled','b');
% line([turnMaxResFreq(6), turnMaxResFreq(6)],[max(max(weight)), min(min(weight))]);
xlabel('Number of turns (Nt)'), ylabel('Induced signal voltage/Coil noise'), zlabel('isol (mm)')
% The freq nr. 18 is about 10kHz
weight = w_temp(:,:,18);
weight2 = w_temp2(:,:,18); 
hSurface2 = surf(N_tM,weight,deltaM);
set(hSurface2,'FaceColor',[0 0 1],'FaceAlpha',0.5);
hSurface = surf(N_tM,weight2,deltaM);
set(hSurface,'FaceColor',[0.7 0.7 1],'FaceAlpha',0.5);
scatter3(N_tM(1,1),weight(1,1),deltaM(1,1),'filled','sk');
scatter3(N_tM(end,1),weight(end,1),deltaM(end,1),'*k');
scatter3(N_tM(choiceInsulation,choiceWire),weight(choiceInsulation,choiceWire),deltaM(choiceInsulation,choiceWire),'filled','b');
% line([turnMaxResFreq(18), turnMaxResFreq(18)],[max(max(weight)), min(min(weight))]);
% The freq nr. 18 is about 21kHz
weight = w_temp(:,:,23);
weight2 = w_temp2(:,:,23); 
hSurface3 = surf(N_tM,weight,deltaM);
set(hSurface3,'FaceColor',[0 1 0],'FaceAlpha',0.5);
hSurface = surf(N_tM,weight2,deltaM);
set(hSurface,'FaceColor',[0.7 1 0.7],'FaceAlpha',0.5);
scatter3(N_tM(1,1),weight(1,1),deltaM(1,1),'filled','sk');
scatter3(N_tM(end,1),weight(end,1),deltaM(end,1),'*k');
scatter3(N_tM(choiceInsulation,choiceWire),weight(choiceInsulation,choiceWire),deltaM(choiceInsulation,choiceWire),'filled','b');
% line([turnMaxResFreq(23), turnMaxResFreq(23)],[max(max(weight)), min(min(weight))]);
% The freq nr. 80 is about 50kHz
weight = w_temp(:,:,29);
weight2 = w_temp2(:,:,29); 
hSurface4 = surf(N_tM,weight,deltaM);
set(hSurface4,'FaceColor',[0 1 1],'FaceAlpha',0.5);
hSurface = surf(N_tM,weight2,deltaM);
set(hSurface,'FaceColor',[0.7 1 1],'FaceAlpha',0.5);
scatter3(N_tM(1,1),weight(1,1),deltaM(1,1),'filled','sk');
scatter3(N_tM(end,1),weight(end,1),deltaM(end,1),'*k');
scatter3(N_tM(choiceInsulation,choiceWire),weight(choiceInsulation,choiceWire),deltaM(choiceInsulation,choiceWire),'filled','b');
% line([turnMaxResFreq(29), turnMaxResFreq(29)],[max(max(weight)), min(min(weight))]);
grid on
title('Pseudo SNR estimation, Body noise not included!')
legend([hSurface1,hSurface2,hSurface3,hSurface4,hLine1_10,hLine1_11],'2kHz','10kHz','21kHz','50kHz','Small diam & isol','Small diam & large isol')