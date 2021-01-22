% Ns = total number of strands.
% NB = number of bunching operations, (agrupations of strands)
% NC = number of cabling,(agrupations of bunches)
function [R_ac_prox,R_ac,Lb,Cb,Rb_dc] = coilRLCcalcLitz_Commercial(N_v,N_l,R_i,rcb,zcb,di_strand,do,dr,dz,freqs,NS,NB,NC)
    % figure(1),plot_wires(rcb,zcb,di,do);                  % Visualize wires for visual confirmation
    l_w = coilLength(rcb);                                  % Calculate the length of the coil
    [Rb_dc] = CoilResistDCLitz_Commercial(l_w,di_strand,NS,NB,NC);  % DC Resistance calculation
    s = do + sqrt(dr^2+dz^2);                               % distance between wire centres
    [R_ac_prox,R_ac] = CoilResistACLitz_Commercial(Rb_dc,rcb,zcb,di_strand,s,freqs,NS,do);    % AC Resistance calculation
    % figure(2),plotResist(Rb_dc,Rb,F,G_2,freqtmp2);        % Visualize Resistance vs. frequency
    [Lb,Lo,M]= CoilInductanceDC(s,rcb,zcb);                 % Inductance (DC)
    di = do + 0.1e-3;
    [Cb]= CoilCapacitanceDC(R_i,di,do,N_v,N_l,dz,dr);     % Capacitance
end