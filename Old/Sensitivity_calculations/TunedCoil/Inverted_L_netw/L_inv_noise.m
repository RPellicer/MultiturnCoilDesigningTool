%% Analyzing noise performance due to an inverted L matching network in between coil and pre-amplifier
function [vn_in, Zth, noisesAtInput, nameNoises] = L_inv_noise(e_n, i_n, e_OutStagNoise, freqs, Rs, Ls, Cs, C1, Qc1, C2, Qc2, G_effective,GainAmp,Ca,Ra,visu)
    % Calculate the thevening impedance seen from the amplifier
%     freqs = logspace(1,6,10000)
    Xl = 1j*2*pi.*freqs.*Ls;
    Xcs = 1./(1j*2*pi.*freqs.*Cs);
    Xc1 = 1./(1j*2*pi.*freqs.*C1);
    Xc2 = 1./(1j*2*pi.*freqs.*C2);
    Xca = 1./(1j*2*pi.*freqs.*Ca);
    %% Analyzing noise from the quality factor of the capacitors (Q = around 100 in the kHz regime) for an inverted L-network
    ESR1 = abs(Xc1)./Qc1;
    ESR2 = abs(Xc2)./Qc2; % <= % Q = abs(Xc)/ESR;
%     Zth = (Rs+Xl+ESR1+Xc1).*(Xc2+ESR2)./(Rs+Xl+ESR1+Xc1+Xc2+ESR2); % This should match with Rmatch
    Zs = (Rs+Xl).*Xcs./(Rs+Xl+Xcs); % In case Cs becomes significant it is taken into consideration
    Zth = (Zs+ESR1+Xc1).*(Xc2+ESR2)./(Zs+Xc1+ESR1+Xc2+ESR2); % This should match with Rmatch
    Za_C2 = 1./(1./(Xc2+ESR2)+1./Xca+1./Ra); % impedance seen after the first capacitor towards the amplifier
    gain_en_C1 = Za_C2./(Za_C2+Zs+Xc1+ESR1);
    % formulas for convertion from series to parallel
    EPR2 = (1+Qc2^2).*ESR2;
    % Noise introduced by non-ideal capacitors in an inverted L network
    en_C1 = en_thermal(ESR1);
    in_C2 = in_thermal(EPR2);
%     en_Rs = en_thermal(Rs);
    en_Rs = en_thermal(Rs);
    % Calculate individually the VOLTAGE noise sources to the input of the amplifier
    ce_ni = abs(e_n).*ones(size(freqs));     % Voltage noise of the amplifier
    ci_ii = abs(i_n.*abs(Zth));          % Current noise of the amplifier
    ce_si = en_Rs.*G_effective; % Thermal noise of the source
    ce_c1i = en_C1.*abs(gain_en_C1);         % Thermal noise of the series matching capacitor (C1)
    ce_c2i = in_C2.*abs(Zth);            % Thermal noise of the parallel tuning capacitor (C2)
    oe_i = e_OutStagNoise.*ones(size(freqs))./GainAmp ;  % This noise dominates with small gains
    ce_tot = sqrt(ce_ni.^2 + ci_ii.^2 + ce_si.^2 + ce_c1i.^2 + ce_c2i.^2 + oe_i.^2); % VOLTAGE noise sources INPUT-REFERRED
    vn_in = ce_tot;
    noisesAtInput = [ce_ni;ci_ii;ce_si;ce_c1i;ce_c2i;oe_i;ce_tot];
    nameNoises = {'ce_ni';'ci_ii';'ce_si';'ce_c1i';'ce_c2i';'oe_i';'ce_tot'};
    %% ################################################## 
    %     VISUALIZATION
    % ##################################################
    if visu == 1
        %% Visualization of the noise seen from the amplifier
        figure(21);
        loglog(freqs,ce_ni,'*')
        hold on
        loglog(freqs,ci_ii,'*')
        loglog(freqs,ce_si,'*')
        loglog(freqs,ce_c1i,'*')
        loglog(freqs,ce_c2i,'*')
        loglog(freqs,oe_i,'*')
        loglog(freqs,ce_tot,'*')
        title('Noise analysis'),xlabel('Freq(Hz)'),ylabel('V/sqrt(Hz)')
        legend('en_a_m_p','in_a_m_p','en_R_s','en_C_1','en_C_2','en_a_m_p_o_u_t_p','en_i_n_-_t_o_t_a_l')
        grid on
        
        % Visualize impedances
        figure(22);
        loglog(freqs,abs(Zth),'r')
        hold on, grid on
        loglog(freqs,real(Zth),'m')
        loglog(freqs,abs(imag(Zth)),'c')
        loglog(freqs,real(Zs),'k')
        loglog(freqs,abs(imag(Zs)),'b')
        loglog(freqs,abs(Za_C2),'g')
        legend('Zth','Zth_r_e_a_l','Zth_i_m_a_g','Zs_r_e_a_l','Zs i_m_a_g','Za_C_2')
    end
end