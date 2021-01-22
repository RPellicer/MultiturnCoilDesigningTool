%% Calculating the capacitors needed for the matching network accounting
% with the ESR (Equivalent Series Resistance) of the capacitor C1
% Camp is the INPUT capacitance of the amplifier. If unknown set Camp = 0;
% [C1, C2, Zth0] = L_inv_estim_cap_ESR(1625, 1e4, 76.4, 55.63e-3, 2e-12, 100)
function [C1, C2, Zth0] = L_inv_estim_cap_ESR(Rmatch, f0, Rs, Ls, Cs, Camp, Qc1, Qc2)
    % Calculations acording to 'lect7_match.pdf' of Berkeley
    w0 = 2*pi.*f0;
    Xls = 1j.*w0.*Ls;
    Xcs = 1/(1j.*w0.*Cs);
    Zs = (Rs+Xls)*Xcs/(Rs+Xls+Xcs);
    Xs_pri = 1j*imag(Zs);
    Rs_pri = real(Zs);
    %% For noise matching
    Rs_ESR = Rs_pri; % Calculate the boosting factor
    for idx = 1:4  % one iteration should be close enough but 4 should give a closer match
        m = Rmatch/Rs_ESR; % Calculate the boosting factor
        % Rpl = RnoiseMatch*Rap/(RnoiseMatch+Rap); % R- Trying to tune for low input impedance amplifiers
        % m = Rpl/Rs; % R- Calculate the boosting factor
        Qreq = sqrt(m-1); % Compute the required circuit Q by (1 + Q2) = m => Q = sqrt(m-1)
        Xs_total = 1j*Qreq*Rs_ESR; % Qreq = Xs_t/Rs;
        Xs_tun1 = Xs_total - Xs_pri;
        Xs_total_p = Xs_total*(1+Qreq^-2); % Equivalent parallel impedance
        Xs_tun2 = -Xs_total_p;

        % Calculate the value of the capacitors needed
        C1 = 1./(1j.*w0*Xs_tun1);
        Ct2_temp = 1./(1j.*w0*Xs_tun2);
        C2 = Ct2_temp-Camp; % Consider the input capacitance of the amplifier

        % Calculate the Thevening impedance seen from the amplifier      
        Xl = 1j*2*pi*f0*Ls;
        Xcs = 1./(1j*2*pi*f0.*Cs);
        Xc1 = 1./(1j*2*pi*f0.*C1);
        Xc2 = 1./(1j*2*pi*f0.*C2);
        %% Analyzing noise from the quality factor of the capacitors (Q = around 100 in the kHz regime) for an inverted L-network
        ESR1 = abs(Xc1)./Qc1; % QvalueManufAdjuster(Q,C,f)
        ESR2 = abs(Xc2)./Qc2; % <= % Q = abs(Xc)/ESR;
    %     Zth = (Rs+Xl+ESR1+Xc1).*(Xc2+ESR2)./(Rs+Xl+ESR1+Xc1+Xc2+ESR2); % This should match with Rmatch
        Zs = (Rs+Xl).*Xcs./(Rs+Xl+Xcs); % In case Cs becomes significant it is taken into consideration
        Zth0 = (Zs+ESR1+Xc1).*(Xc2+ESR2)./(Zs+Xc1+ESR1+Xc2+ESR2); % This should match with Rmatch
        Rs_ESR = Rs_pri+1/(2*pi*f0*C1*Qc1);          
    end 
%     if real(C1)<0 || real(C2)<0 || imag(Zs)<0 || Rs_ESR >= Rmatch % otherwise it sees a capacitance and can't be tuned
%         C1 = 1e-17;
%         C2 = 1e-17;
%         C1 = NaN;
%         C2 = NaN;
%     end
%     if real(C1)<0 % This happens when required Q is higher than the one of the coil
%         C1= 1e18; % This means that the capacitor is a short cirquit. The ESR of this capacitor is very small => Rc_ESR = 1/(w*C*Q)
%     end
%     if real(C2)<0
%         C2= 1e-18; % This means that the capacitor is an open cirquit.
%     end
end