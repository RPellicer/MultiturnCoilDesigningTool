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

function [sensScore,probe] = sensitivityI2V(tProbe,visu)
    e_n = tProbe.Amp.e_n_glob;
    i_n = tProbe.Amp.i_n_glob; 
    e_OutStagNoise = tProbe.Amp.e_OutStagNoise_glob;
%     Camp = tProbe.Amp.Camp_glob;
%     Ramp = tProbe.Amp.Ramp_glob;
%     GainAmp = tProbe.Amp.GainAmp_glob;
    Rf = tProbe.Amp.Rf_glob;
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
    % possition of the wires
    [rcb,zcb] = XY_wire_ort(N_v,N_l,R_i,do,dr,dz);        % Get the location of each wire
    % Resistance, inductance and capacitance of the coil
    if (tProbe.Coil.LitzStrandNr_glob == 1)
        [Rs,Ls,Cs,Rs_dc] = coilRLCcalc(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,Freqs);
    else
        [Rs,Ls,Cs,Rs_dc] = coilRLCcalcLitz(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,Freqs,NS,NB,NC);
    end
    l_w = coilLength(rcb);
    Rs0 = Rs(ceil(length(Rs)/2));
    % Voltage noise at the input
    [vn_in, noisesAtInput, nameNoises] = ampNoiseI2V(i_n,e_n,e_OutStagNoise,Rf,Rs,Ls,Cs,Freqs,visu); % V/sqrt(Hz)
    % Field transfer ratio (emf(V)/T)
    if NearField
        emf_per_T = emfNear(Freqs,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp);
    else
        emf_per_T = emfFar(Freqs,rcb); % This is normalized to VOLTS PER TESLA (V/T)
    end
    % for the Tina program (uses an inductor to fake the frequency dependence)
    % emf_per_T=1j*2*pi*freqs*L => L = emf_per_T/(1j*2*pi*freqs)
    % Convert the input voltage noise to magnetic field noise to know the
    % sensitivity to the field
    sens = vn_in./abs(emf_per_T);  
    sensPerNoise = noisesAtInput./repmat(emf_per_T,size(noisesAtInput,1),1);
%     sensMean = trapz(Freqs,sens)/(Freqs(end)-Freqs(1));
    sensMean = mean(sens); % The noise floor is averaged
    sens_peak = min(sens);

%% VISUALIZATION
    if visu == 1
    % Visualization of the sensitivity
        figure(10);
        loglog(Freqs,sens,'-*')
        title('sensitivity of the coil+amplifier'),xlabel('Freq(Hz)'),ylabel('Field sensitivity T/sqrt(Hz)')
        grid on, axis tight
%     % Visualization of the vn_in
%         figure(11);
%         loglog(freqs,vn_in,'-*')
%         title('Noise of the coil+amplifier'),xlabel('Freq(Hz)'),ylabel('Total electronic noise V/sqrt(Hz)')
%         grid on, axis tight
%     % Visualization of the emf at the amplifier
%         figure(12);
%         loglog(freqs,emf_at_amp,'-*')
%         title('emf'),xlabel('Freq(Hz)'),ylabel('emf at the amplifier V/sqrt(Hz)')
%         grid on, axis tight
%     % Visualization of the Zth seen from the amplifier
%         figure(13);
%         loglog(freqs,abs(Zth),'-*')
%         title('Zth'),xlabel('Freq(Hz)'),ylabel('ohms')
%         grid on, axis tight
     % Plot the wires (cross section)
        figure(14);
        plot_wires(rcb,zcb,di,do)
        title('Side cut view of the coil profile')
    end
    %%     Save new variables of interest into structures 
    % Adding the inputed data to 'probe'
    probe.Amp = tProbe.Amp;
    probe.Coil = tProbe.Coil;
    probe.Freqs = tProbe.Freqs;
    % Adding the results
    % Avoid coil designs that are too close to the self-resonance
    % frequency
    probe.Coil.SelfRes = 1/(2*pi*sqrt(Ls*Cs));
    distFreq = 10;
    if probe.Coil.SelfRes > f0*distFreq;
        probe.Sens.PenalyzeScaling = tProbe.Coil.PenalyzeScaling;
    else
%         probe.Sens.PenalyzeScaling = tProbe.Coil.PenalyzeScaling;
        probe.Sens.PenalyzeScaling = abs(1+((f0*distFreq-probe.Coil.SelfRes)*1/(f0*distFreq))) + tProbe.Coil.PenalyzeScaling;
    end
    sensScore = sensMean*probe.Sens.PenalyzeScaling;
    probe.Sens.SensScore = sensScore;
    probe.Sens.SensMean = sensMean;
    probe.Sens.SensRange = sens;
    probe.Sens.Sens_peak =sens_peak;
    probe.Sens.SensPerNoise = sensPerNoise;
    probe.Sens.NameNoises = nameNoises;
    probe.Freqs.Freqs = Freqs;
    probe.Coil.Rs0 = Rs0;
    probe.Coil.Rs = Rs;
    probe.Coil.Rs_dc = Rs_dc;
    probe.Coil.Ls = Ls;
    probe.Coil.Cs = Cs;
    probe.Coil.Emf_per_T = emf_per_T;
    probe.Coil.PspiceFreqDep = abs(mean(probe.Coil.Emf_per_T./(2*pi*probe.Freqs.Freqs)));
    probe.Coil.WireLength = l_w;
    probe.Type = 'Non-tuned';
end