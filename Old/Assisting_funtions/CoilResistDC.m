


function R=CoilResistDC(l,di)
% DC Resistance calculation. Everything in standard international units
% ro_c = 1.7e-8; % resistivity of the wire
R = (4*1.7e-8/(di^2))*(l/pi); % Total DC resistance of the cable
end