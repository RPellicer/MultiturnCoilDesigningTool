% Ns = total number of strands.
% NB = number of bunching operations, (agrupations of strands)
% NC = number of cabling,(agrupations of bunches)
function [Rb,Lb,Cb,Rb_dc] = coilRLCcalcLitz(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,freqs,NS,NB,NC)
    % figure(1),plot_wires(rcb,zcb,di,do);                  % Visualize wires for visual confirmation
    l_w = coilLength(rcb);                                  % Calculate the length of the coil
    [Rb_dc,di_strand] = CoilResistDCLitz(l_w,di,NS,NB,NC);  % DC Resistance calculation
    s = do + sqrt(dr^2+dz^2);                               % distance between wire centres
    [Rb,F,G_2] = CoilResistACLitz(Rb_dc,rcb,zcb,di_strand,s,freqs,NS);    % AC Resistance calculation
    % figure(2),plotResist(Rb_dc,Rb,F,G_2,freqtmp2);        % Visualize Resistance vs. frequency
    [Lb,Lo,M]= CoilInductanceDC(s,rcb,zcb);                 % Inductance (DC)
    [Cb]= CoilCapacitanceDC(R_i,di,do,N_v,N_l,dz,dr);     % Capacitance
end