function [Cs]= CoilCapacitanceDC(Ri,di,do,Nv,Nl,Ksi_z,Ksi_r)
    eo = 8.854187817e-12; % Vacuum electrical permittivity F/m (farads per metre)
    er = 4.4; % Martinez used 2.7;
    % di = Diameter of the wire without coating in meters
    % do = Diameter of the wire with coating in meters
    dc = Nv*do+(Nv-1)*Ksi_z;% Coil thickness
    % From Martinez et al, 2014 (eq. nr. 44)
    d_ = 1.26*do-1.15*di+Ksi_z;
    % From Martinez et al, 2014 (eq. nr. 46)
%     Cs =(8*pi*eo*er*dc*(Nl-1)/(6*d_*Nl^2))*(2*Ri+do-Nl*Ksi_r);
    if Nl == 1
        e0 = 8.85418e-12; % F/m
        Cs  = (pi^2)*Ri*2*e0/log(do/di+sqrt((do/di)^2-1)); % From Grandi-1999 eq. 12
%         t = do-di;
%         F = (do/di)/(1+t/(di/2))^(1-1/1); % From Grandi-1999 eq. 13
%         Cs2 = (pi^2)*Ri*2*e0/log(F+sqrt(F^2-(1+t/(di/2))^(2/1))); % From Grandi-1999 eq. 14        
    else
        Cs =(8*pi*eo*er*dc*(Nl-1)/(6*d_*Nl^2))*(2*Ri+do-Nl*Ksi_r);
    end
end

% refs:
% Martinez-2014-On Evaluation of Inductance, DC Resistance, and Capacitance of Coaxial Inductors at Low Frequencies
% Grandi-1999-Stray capacitances of single-layer solenoid air-core inductors