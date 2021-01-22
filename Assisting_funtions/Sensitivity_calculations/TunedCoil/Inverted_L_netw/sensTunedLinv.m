%% This function calculates the sensitivity (T/sqrt(Hz)) of a surface coil with rectangular cross-section

% di [METERS] cooper outer diameter
% do [METERS] wire outer diameter insulation included
% N_v  Number of turns per layer
% N_l Number of layers
% dr axial filling factor
% dz azimut filling factor
% R_i Coil internal radius

% Asume constant values for optimization: i_n,e_n,e_OutStagNoise,Rf,dr,dz,freqs
% Optimize for: N_v,N_l,R_i,di,do
% ONLY ONE FREQUENCY
% sensTunedLinv(1.3e-9,0.8e-12,90e-9,42,35,10e-3,0.2e-3,0.23e-3,0,0,1e4,0,2e-12,60e6,100,1000)
% Show the whole BW
% [sensMean,sens,sensf0,freqs] = sensTunedLinv(1.3e-9,0.8e-12,90e-9,42,35,10e-3,0.2e-3,0.23e-3,0,0,3.3e3,6e3,2e-12,60e6,100,100,1000);
% [sensMean,sens,sens_peak,freqs,C1,C2,Rs,Ls,Cs,emf_per_T,C1_ESR, C2_ESR] = sensTunedLinv(e_n,i_n,e_OutStagNoise,N_v,N_l,R_i,di,do,dr,dz,f0,bw,Camp,Ramp,Qcap1ESR,Qcap2ESR,GainAmp,visu)
% 
function [sensScore, probe] = sensTunedLinv(tProbe,visu)
%%     Structure variables are deployed for visivility
    e_n = tProbe.Amp.e_n_glob;
    i_n = tProbe.Amp.i_n_glob; 
    e_OutStagNoise = tProbe.Amp.e_OutStagNoise_glob;
    Camp = tProbe.Amp.Camp_glob;
    Ramp = tProbe.Amp.Ramp_glob;
    GainAmp = tProbe.Amp.GainAmp_glob;
    N_v = tProbe.Coil.N_v;
    N_l = tProbe.Coil.N_l;
    R_i = tProbe.Coil.Rin;
    di = tProbe.Coil.di;
    do = tProbe.Coil.do;
    dr = tProbe.Coil.dr;
    dz = tProbe.Coil.dz;
    
    NS = tProbe.Coil.NS; % total number of strands.
    NB = tProbe.Coil.NB; % number of bunching operations, (agrupations of strands)
    NC = tProbe.Coil.NC; % number of cabling,(agrupations of bunches)
    bw = tProbe.Freqs.bw_glob;
    f0_Offset = tProbe.Freqs.f0_Offset;
    f0_cent = tProbe.Freqs.f0_glob;
    f0 = f0_cent + f0_Offset;
    Qcap1ESR = tProbe.MatchNetw.Qcap1ESR_glob;
    Qcap2ESR = tProbe.MatchNetw.Qcap2ESR_glob;
    NearField = tProbe.Coil.emfNear;
    Bp =tProbe.Coil.Bp;
    r_ROI = tProbe.Coil.r_ROI;
    hBase_ROI = tProbe.Coil.hBase_ROI;
    hfinal_ROI = tProbe.Coil.hfinal_ROI; 

    if bw == 0
        Freqs = f0;
    else
        Freqs = linspace(f0_cent-bw/2,f0_cent+bw/2,99);
    end
    %%    
    % possition of the wires
    [rcb,zcb] = XY_wire_ort(N_v,N_l,R_i,do,dr,dz);        % Get the location of each wire
    % Resistance, inductance and capacitance of the coil
    if (tProbe.Coil.LitzStrandNr_glob == 1)
        [Rs,Ls,Cs,Rs_dc] = coilRLCcalc(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,Freqs);
    else
        [Rs,Ls,Cs,Rs_dc] = coilRLCcalcLitz(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,Freqs,NS,NB,NC);
    end
    Rs0 = Rs(ceil(length(Rs)/2));
    Cs = Cs + 1e-15; % Adding a minimum of capacitance avoids errors due to Cs = 0
    % Calculate the matching network to noise match the amplifier
    kb = 1.38e-23; %k is Boltzmann’s constant (J/K)
    T = 300; % room temperature in kelvin
    % Search the best matching: It should be in between e_rs =e_i and e_n = e_i; 
%     RmatchMax = (4*kb*T/(i_n^2)); % in*Rs = e_rs => in^2*Rs^2 = 4*k*T*Rs +> Rs = 4*k*T/in^2; Noise from the resistor and current noise are equal
    RmatchOfficial = e_n/i_n;
    RmatchMax = RmatchOfficial * 100;
    RmatchMin = RmatchOfficial / 100;  % e_n/i_n is the Usual noise match of the amplifier
    RmatchSearchPointNr = 21; % the theoretical optimal e_n/i_n is in the 11th possition
    Rmatch = logspace(log10(RmatchMin),log10(RmatchMax),RmatchSearchPointNr); % Desired virtual resistance exposed to the amplifier
    Rmatch(11) = RmatchOfficial; % Sometimes the number is not precise
    % Estimate (rougly) the vaues of the CAPACITORS
    for k = 1:length(Rmatch)
        [C1, C2, ~] = L_inv_estim_cap_ESR(Rmatch(k), f0, Rs0, Ls, Cs, Camp, Qcap1ESR, Qcap2ESR);
        [tempGain, ~, ~] = freqRespLinv_ESR(Rs,Ls,Cs,C1,C2,Camp,Ramp,Freqs,Qcap1ESR,Qcap2ESR);
        G_effective = abs(tempGain);
        % Estimate the noise at the input of the amplifier
        [vn_in, ~, ~, ~] = L_inv_noise(e_n, i_n, e_OutStagNoise, Freqs, Rs, Ls, Cs, C1, Qcap1ESR, C2, Qcap2ESR, G_effective, GainAmp, Camp, Ramp, 0); % V/sqrt(Hz)   
        % Field transfer ratio (emf(V)/T)
        if NearField
            emf_per_T = emfNear(Freqs,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp);
        else
            emf_per_T = emfFar(Freqs,rcb); % This is normalized to VOLTS PER TESLA (V/T)
        end            %     TinaFreqDep = emf_per_T/(2*pi*freqs);
        emf_at_amp = abs(emf_per_T).*G_effective; % normalized to VOLTS PER TESLA (V/T)
        % Convert the input voltage noise to magnetic field noise to know the
        % sensitivity to the field
        sens = vn_in./emf_at_amp;
        tSensMean(k) = mean(sens); % The noise floor is averaged
%         tSensMean(k) = nansum(sens)/length(sens); % The noise floor is averaged
%             if (C1==1e+18 || C2==1e+18) % Added to prevent using a single capacitor tuning. This is not a proble other than for very low frequencies where coils would get large amount of loops.
%                 tSensMean(k) = 1e+18;
%             end
        if ((imag(C1)~=0)||(imag(C2)~=0)||(C1>1e-3)||(C2>1e-3)) % non realistic capacitor values
            tSensMean(k) = 1e+18;
        end
    end
    [~,ind] = min(tSensMean);
    if visu == 1
        figure(9);
        loglog(Rmatch,tSensMean,'r','Marker','*');
        hold on, grid on, xlabel('Rmatch'), ylabel('Sensitivity'), title('Optimal Rmatch')
        loglog(RmatchOfficial,tSensMean(find(Rmatch == RmatchOfficial)),'b','Marker','*')
    end
%     else
%         Rmatch = logspace(log10(RmatchMin),log10(RmatchMax),90); % Desired virtual resistance exposed to the amplifier
%         % Estimate (rougly) the vaues of the CAPACITORS
%         for k = 1:length(Rmatch)
%             [C1, C2, ~] = L_inv_estim_cap_ESR(Rmatch(k), f0, Rs0, Ls, Cs, Camp, Qcap1ESR, Qcap2ESR);
%             [tempGain, ~, ~] = freqRespLinv_ESR(Rs,Ls,Cs,C1,C2,Camp,Ramp,Freqs,Qcap1ESR,Qcap2ESR);
%             G_effective = abs(tempGain);
%             % Estimate the noise at the input of the amplifier
%             [vn_in, ~, ~, ~] = L_inv_noise(e_n, i_n, e_OutStagNoise, Freqs, Rs, Ls, Cs, C1, Qcap1ESR, C2, Qcap2ESR, G_effective, GainAmp, Camp, Ramp, 0); % V/sqrt(Hz)   
%             % Field transfer ratio (emf(V)/T)
%             if NearField
%                 emf_per_T = emfNear(Freqs,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp);
%             else
%                 emf_per_T = emfFar(Freqs,rcb); % This is normalized to VOLTS PER TESLA (V/T)
%             end            %     TinaFreqDep = emf_per_T/(2*pi*freqs);
%             emf_at_amp = abs(emf_per_T).*G_effective; % normalized to VOLTS PER TESLA (V/T)
%             % Convert the input voltage noise to magnetic field noise to know the
%             % sensitivity to the field
%             sens = vn_in./emf_at_amp;
%             tSensMean(k) = nansum(sens)/length(sens); % The noise floor is averaged
%             if (C1==1e+18 || C2==1e+18) % Added to prevent using a single capacitor tuning. This is not a problem other than for very low frequencies where coils would get large amount of loops.
%                 tSensMean(k) = 1e+18;
%             end
%         end
%         [~,ind] = min(tSensMean);
%         RmatchT = logspace(log10(RmatchMin),log10(RmatchMax),RmatchSearchPointNr); % Desired virtual resistance exposed to the amplifier
%         % Estimate (rougly) the values of the CAPACITORS
%         for k = 1:length(RmatchT)
%             [C1, C2, ~] = L_inv_estim_cap_ESR(RmatchT(k), f0, Rs0, Ls, Cs, Camp, Qcap1ESR, Qcap2ESR);
%             [tempGain, ~, ~] = freqRespLinv_ESR(Rs,Ls,Cs,C1,C2,Camp,Ramp,Freqs,Qcap1ESR,Qcap2ESR);
%             G_effective = abs(tempGain);
%             % Estimate the noise at the input of the amplifier
%             [vn_in, ~, ~, ~] = L_inv_noise(e_n, i_n, e_OutStagNoise, Freqs, Rs, Ls, Cs, C1, Qcap1ESR, C2, Qcap2ESR, G_effective, GainAmp, Camp, Ramp, 0); % V/sqrt(Hz)   
%             % Field transfer ratio (emf(V)/T)
%             if NearField
%                 emf_per_T = emfNear(Freqs,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp);
%             else
%                 emf_per_T = emfFar(Freqs,rcb); % This is normalized to VOLTS PER TESLA (V/T)
%             end
%             emf_at_amp = abs(emf_per_T).*G_effective; % normalized to VOLTS PER TESLA (V/T)
%                 %% Visualization of the emf at the amplifier
% %             figure(12);
% %             loglog(Freqs,emf_at_amp,'-*')
% %             hold on
% %             title('emf'),xlabel('Freq(Hz)'),ylabel('emf at the amplifier V/sqrt(Hz)')
%             % Convert the input voltage noise to magnetic field noise to know the
%             % sensitivity to the field
%             sensT = vn_in./emf_at_amp;
%             tSensMeanT(k) = nansum(sensT)/length(sensT); % The noise floor is averaged
%         end
% %         [~,ind] = min(tSensMeanT);
%         figure(9);
%         loglog(Rmatch,tSensMean);
%         hold on, grid on, xlabel('Rmatch'), ylabel('Sensitivity'), title('Optimal Rmatch')
%         loglog(RmatchT,tSensMeanT,'r*')
%     end
    [C1, C2, ~] = L_inv_estim_cap_ESR(Rmatch(ind), f0, Rs0, Ls, Cs, Camp, Qcap1ESR, Qcap2ESR);
% %     This manipulation untuned the coil if the values of the capacitors
    [tempGain, C1_ESR, C2_ESR] = freqRespLinv_ESR(Rs,Ls,Cs,C1,C2,Camp,Ramp,Freqs,Qcap1ESR,Qcap2ESR);
    G_effective = abs(tempGain);
    % Estimate the noise at the input of the amplifier
    [vn_in, ~, noisesAtInput,nameNoises] = L_inv_noise(e_n, i_n, e_OutStagNoise, Freqs, Rs, Ls, Cs, C1, Qcap1ESR, C2, Qcap2ESR, G_effective, GainAmp, Camp, Ramp, visu); % V/sqrt(Hz)   
    % Field transfer ratio (emf(V)/T)
    if NearField
        emf_per_T = emfNear(Freqs,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp);
    else
        emf_per_T = emfFar(Freqs,rcb); % This is normalized to VOLTS PER TESLA (V/T)
    end    %     TinaFreqDep = emf_per_T/(2*pi*freqs);
    emf_at_amp = abs(emf_per_T).*G_effective; % normalized to VOLTS PER TESLA (V/T)
    % Convert the input voltage noise to magnetic field noise to know the
    % sensitivity to the field
    sens = vn_in./emf_at_amp;
    sensPerNoise = noisesAtInput./repmat(emf_at_amp,size(noisesAtInput,1),1);
    sensMean = nansum(sens)/length(sens); % The noise floor is averaged
%     sensMean = trapz(Freqs,sens)/(Freqs(end)-Freqs(1));
    sens_peak = min(sens);
    %% ################################################## 
    %     VISUALIZATION
    % ##################################################
    Freqs2 = linspace(f0_cent-bw-100,f0_cent+bw+100,199); % For visualising porposes only
    if (tProbe.Coil.LitzStrandNr_glob == 1)
        [Rs2,~,~,~] = coilRLCcalc(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,Freqs2);
    else
        [Rs2,~,~,~] = coilRLCcalcLitz(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,Freqs2,NS,NB,NC);
    end
    l_w = coilLength(rcb);
    [tempGain2, C1_ESR, C2_ESR] = freqRespLinv_ESR(Rs2,Ls,Cs,C1,C2,Camp,Ramp,Freqs2,Qcap1ESR,Qcap2ESR);
    G_effective2 = abs(tempGain2);
    % Estimate the noise at the input of the amplifier
    [vn_in2, ~, noisesAtInput2,nameNoises] = L_inv_noise(e_n, i_n, e_OutStagNoise, Freqs2, Rs2, Ls, Cs, C1, Qcap1ESR, C2, Qcap2ESR, G_effective2, GainAmp, Camp, Ramp, visu); % V/sqrt(Hz)   
    % Field transfer ratio (emf(V)/T)
    if NearField
        emf_per_T2 = emfNear(Freqs2,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp);
    else
        emf_per_T2 = emfFar(Freqs2,rcb); % This is normalized to VOLTS PER TESLA (V/T)
    end    %     TinaFreqDep = emf_per_T/(2*pi*freqs);
    emf_at_amp2 = abs(emf_per_T2).*G_effective2; % normalized to VOLTS PER TESLA (V/T)
    % Convert the input voltage noise to magnetic field noise to know the
    % sensitivity to the field
    sens2 = vn_in2./emf_at_amp2;
%     sensPerNoise2 = noisesAtInput2./repmat(emf_at_amp2,size(noisesAtInput2,1),1);

    if visu == 1
    %% Visualization of the sensitivity
        figure(10);
        loglog(Freqs2,sens2,'-*')
        hold on
        loglog(Freqs,sens,'-*')
        title('sensitivity of the coil+amplifier'),xlabel('Freq(Hz)'),ylabel('Field sensitivity T/sqrt(Hz)')
        grid on, axis tight
    %     %% Visualization of the vn_in
    %         figure(11);
    %         loglog(freqs,vn_in,'-*')
    %         title('Noise of the coil+amplifier'),xlabel('Freq(Hz)'),ylabel('Total electronic noise V/sqrt(Hz)')
    %         grid on, axis tight
    %% Visualization of the emf at the amplifier
        figure(12);
        loglog(Freqs2,emf_at_amp2,'-*')
        hold on
        loglog(Freqs,emf_at_amp,'-*')
        title('emf'),xlabel('Freq(Hz)'),ylabel('emf at the amplifier V/sqrt(Hz)')
        grid on, axis tight
    %     %% Visualization of the Zth seen from the amplifier
    %         figure(13);
    %         loglog(freqs,abs(Zth),'-*')
    %         title('Zth'),xlabel('Freq(Hz)'),ylabel('ohms')
    %         grid on, axis tight
     %% Plot the wires (cross section)
        figure(14);
        plot_wires(rcb,zcb,di,do)
        if NearField
            line([hBase_ROI,hBase_ROI,hfinal_ROI,hfinal_ROI],[0,r_ROI,r_ROI,0]);%             r_ROI,hBase_ROI,hfinal_ROI
        end
        title('Side cut view of the coil profile')
    end
    %%     Save new variables of interest into structures 
    % Adding the inputed data to 'probe'
    probe.Amp = tProbe.Amp;
    probe.Coil = tProbe.Coil;
    probe.Freqs = tProbe.Freqs;
    probe.MatchNetw = tProbe.MatchNetw;
    % Adding the results
    % Avoid coil designs that are too close to the self-resonance
    % frequency
    probe.Coil.SelfRes = 1/(2*pi*sqrt(Ls*Cs));
    distFreq = 10;
    if probe.Coil.SelfRes > f0*distFreq
        probe.Sens.PenalyzeScaling = tProbe.Coil.PenalyzeScaling;
    else
%         probe.Sens.PenalyzeScaling = tProbe.Coil.PenalyzeScaling;
        probe.Sens.PenalyzeScaling = abs(1+((f0*distFreq-probe.Coil.SelfRes)*1/(f0*distFreq))) + tProbe.Coil.PenalyzeScaling;
    end
%     avoid designs with non-realistic capacitor values
    if ((imag(C1)~=0)||(imag(C2)~=0)||(C1>1e-3))
        probe.Sens.PenalyzeScaling = probe.Sens.PenalyzeScaling * 100;
    end
    sensScore = sensMean*probe.Sens.PenalyzeScaling;
    probe.Sens.SensScore = sensScore;
    probe.Sens.SensMean = sensMean;
    probe.Sens.SensRange = sens;
    probe.Sens.Sens_peak = sens_peak;
    probe.Sens.RmatchSens(1,:) = tSensMean;
    probe.Sens.RmatchSens(2,:) = Rmatch;
    probe.Sens.SensPerNoise = sensPerNoise;
    probe.Sens.NameNoises = nameNoises;
    probe.Freqs = tProbe.Freqs;
    probe.Freqs.Freqs = Freqs;
    probe.Coil.Rs0 = Rs0;
    probe.Coil.Rs = Rs;
    probe.Coil.Rs_dc = Rs_dc;
    probe.Coil.Ls = Ls;
    probe.Coil.Cs = Cs;
    probe.Coil.Emf_per_T = emf_per_T;
    probe.Coil.WireLength = l_w;
    probe.MatchNetw.C1 = C1;
    probe.MatchNetw.C2 = C2;
    probe.MatchNetw.C1_ESR = C1_ESR;
    probe.MatchNetw.C2_ESR = C2_ESR;   
    probe.MatchNetw.Rmatch = Rmatch(ind);
    probe.MatchNetw.RmatchEval = Rmatch;
    probe.Coil.PspiceFreqDep = abs(mean(probe.Coil.Emf_per_T./(2*pi*probe.Freqs.Freqs)));
    probe.Type = 'Tuned';
    if (probe.Coil.Rin ~= tProbe.Coil.Rin_glob) && NearField
       warning('1 sensTuned: inner radius od the resultant coil (Rin) is different to specified inner radius (Rin_glob). Please, verify it complies your specifications')
    end
end