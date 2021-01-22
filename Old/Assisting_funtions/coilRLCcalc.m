function [Rb,Lb,Cb,Rb_dc] = coilRLCcalc(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,freqs)
    % figure(1),plot_wires(rcb,zcb,di,do);                  % Visualize wires for visual confirmation
    l_w = coilLength(rcb);                                  % Calculate the length of the coil
    Rb_dc = CoilResistDC(l_w,di);                           % DC Resistance calculation
    s = do + sqrt(dr^2+dz^2);                               % distance between wire centres
    [Rb,F,G_2] = CoilResistAC(Rb_dc,rcb,zcb,di,s,freqs);    % AC Resistance calculation
    % figure(2),plotResist(Rb_dc,Rb,F,G_2,freqtmp2);        % Visualize Resistance vs. frequency
    [Lb,Lo,M]= CoilInductanceDC(s,rcb,zcb);                 % Inductance (DC)
    [Cb]= CoilCapacitanceDC(R_i,di,do,N_v,N_l,dz,dr);     % Capacitance
end