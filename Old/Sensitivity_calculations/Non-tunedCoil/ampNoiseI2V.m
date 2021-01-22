% Estimate Noise at the entrance of the Non-tuned TRANS MPEDANCE AMPLIFIER
function [vn_in,noisesAtInput,nameNoises] = ampNoiseI2V(i_n,e_n,e_OutStagNoise,Rf,Rcoil,Lcoil,Ccoil,freqs,visu)
    %% Add noise model here for the TRANS-IMPEDANCE AMPLIFIER
    % Amplifier speCcoil. The noise model is taken from "Calculating noise figure in op amps" (Texas Instrument) which matches with 
    % % *The circuits and filteRcoil handbook* FIGURE 15.27 in Chapter 15
    % i_n = 0.8e-12; %INA217
    % e_n = 1.2e-9;  %INA217
    % e_OutStagNoise = 90e-9; %INA217
    w = 2*pi*freqs;
%     Zsen = (Rcoil + 1j*w.*Lcoil)./((1-Lcoil*Ccoil*w.^2)+1j*w.*Rcoil*Ccoil); % Z_Thevenin Acording to Timofeeva et al, 2011 -"Differential Search CoiLcoil Based MagnetometeRcoil: Conditioning, Magnetic Sensitivity, Spatial Resolution"
    Zs = Rcoil + 1j*w.*Lcoil; % Impedance of the branch (inductance + resistance). The other branch (i.e., the inductance) doesn't account in this
    Zf = Rf;
    
    % Define amplifier gains
    sigGain = abs(Zf./Zs);
    noiseGain = abs(1+(Zf./Zs));
    
    % % Calculate individually the VOLTAGE noise sources INPUT-REFERRED
    ce_ni = abs(e_n).*ones(size(freqs));     % Voltage noise of the amplifier
    ci_ii = abs(i_n.*Zs);                    % Current noise of the amplifier
    ce_si = en_thermal(Rcoil);                % Thermal noise of the source
%     ce_fi = en_thermal(Rcoil)/sigGain;      % Thermal noise of the feedback impedance
    oe_i = e_OutStagNoise.*ones(size(freqs))./sigGain ;  % This noise dominates with small gains
%     ce_in = sqrt(ce_ni^2 + ci_ii.^2 + ce_si.^2 + ce_fi.^2 + oe_i.^2); % VOLTAGE noise sources INPUT-REFERRED
    ce_tot = sqrt(ce_ni.^2 + ci_ii.^2 + ce_si.^2 + oe_i.^2); % VOLTAGE noise sources INPUT-REFERRED
    
    % The noise can be calculated to the output by multiplying the noise
    % refered to the input as in "Calculating noise figure in op amps".
    vn_in = ce_tot;
    noisesAtInput = [ce_ni; ci_ii; ce_si; oe_i; ce_tot];
    nameNoises = {'ce_ni'; 'ci_ii'; 'ce_si'; 'oe_i'; 'ce_tot'};

    if visu == 1
        %% Visualization of the noise seen from the amplifier
        figure(21);
        loglog(freqs,ce_ni,'*')
        hold on
        loglog(freqs,ci_ii,'*')
        loglog(freqs,ce_si,'*')
        loglog(freqs,oe_i,'*')
        loglog(freqs,ce_tot,'*')
        title('Noise analysis'),xlabel('Freq(Hz)'),ylabel('V/sqrt(Hz)')
        legend('en_a_m_p','in_a_m_p','en_R_s','en_a_m_p_o_u_t_p','en_i_n_-_t_o_t_a_l')
        grid on
%         % Visualize impedances
%         figure(22);
%         loglog(freqs,Rcoil,'k')
%         hold on, grid on
%         loglog(freqs,abs(1j*w.*Lcoil),'b')
%         legend('Zs_real','Zs imag')
    end
end
