% Rough estimation of the capacitors needed for the matching network
% Camp is the INPUT capacitance of the amplifier. If unknown set Camp = 0;
function [C1, C2, Zth0] = L_inv_estim_cap(Rmatch, f0, Rs, Ls, Camp)
    % Calculations acording to 'lect7_match.pdf' of Berkeley
    w0 = 2*pi.*f0;
    Xs = 1j.*w0.*Ls;
    %% For noise matching
    m = Rmatch/Rs; % Calculate the boosting factor
    % Rpl = RnoiseMatch*Rap/(RnoiseMatch+Rap); % R- Trying to tune for low input impedance amplifiers
    % m = Rpl/Rs; % R- Calculate the boosting factor
    Qreq = sqrt(m-1); % Compute the required circuit Q by (1 + Q2) = m => Q = sqrt(m-1)
    Xs_total = 1j*Qreq*Rs; % Qreq = Xs_t/Rs;
    Xs_tun1 = Xs_total - Xs;
    Xs_total_p = Xs_total*(1+Qreq^-2); % Equivalent parallel impedance
    Xs_tun2 = -Xs_total_p;

    % Calculate the value of the capacitors needed
    C1 = 1./(1j.*w0*Xs_tun1);
    Ct2_temp = 1./(1j.*w0*Xs_tun2);
    C2 = Ct2_temp-Camp; % Consider the input capacitance of the amplifier

    % Calculate the Thevening impedance seen from the amplifier
    Xl = 1j*2*pi*f0*Ls;
    Xc1 = 1./(1j*2*pi*f0.*C1);
    Xc2 = 1./(1j*2*pi*f0.*C2);
    Zth0 = (Rs+Xl+Xc1).*Xc2./(Rs+Xl+Xc1+Xc2); % This should match with Rmatch
end