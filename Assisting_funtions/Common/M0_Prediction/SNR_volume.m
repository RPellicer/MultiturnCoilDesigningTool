%% This function calculates the SNR (Vsignal/Vnoise) of a surface coil with 
% rectangular cross-section and a cylindrical and coaxial sample volume

% di [METERS] cooper outer diameter
% do [METERS] wire outer diameter insulation included
% N_v  Number of turns per layer
% N_l Number of layers
% dr axial filling factor
% dz azimut filling factor
% R_i Coil internal radius

% Asume constant values for optimization: i_n,e_n,e_OutStagNoise,Rf,dr,dz,freqs
% Optimize for: N_v,N_l,R_i,di,do

function [SNR,vn_in,emf] = SNR_volume(i_n,e_n,e_OutStagNoise,Rf,N_v,N_l,R_i,di,do,dr,dz,freqs,r_ROI,hBase_ROI,hfinal_ROI,Bp)
    % possition of the wires
    [rcb,zcb] = XY_wire_ort(N_v,N_l,R_i,do,dr,dz);        % Get the location of each wire
%     % Plot the wires (cross section)
%     plot_wires(rcb,zcb,di,do)
    % Resistance, inductance and capacitance of the coil
    [Rb,Lb,Cb] = coilRLCcalc(N_v,N_l,R_i,rcb,zcb,di,do,dr,dz,freqs);
    % Voltage noise at the input
    vn_in = coilAmpINANoiseTI(i_n,e_n,e_OutStagNoise,Rf,Rb,Lb,Cb,freqs);
    % Field transfer ratio (emf(V)/T)
%     emf = emfFar(freqs,rcb);
    emf = emfNear(freqs,rcb,zcb,r_ROI,hBase_ROI,hfinal_ROI,Bp); % This is the emf (V) generated from the specific ROI.
    % Convert the input voltage noise to magnetic field noise to know the
    % sensitivity to the field
    SNR = abs(emf)./vn_in;
%     sensTot = mean(sens);
end