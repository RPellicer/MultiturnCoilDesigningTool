% NB = number of bunching operations, (agrupations of strands)
% NC = number of cabling,(agrupations of bunches)
% Ns = total number of strands.
% TO BE USED WITH strand number above 100
function [R_DC]=CoilResistDCLitz_Commercial(l,di_strand,NS,NB,NC)
% DC Resistance calculation. Everything in standard international units
% Also see "https://www.newenglandwire.com/products/litz-wire-and-formed-cables/theory"
% ro_c = 1.7e-8; % resistivity of the wire
% FORMULA TAKEN FROM: http://www.elektrisola.com/hf-litz-wire/terminology-basics/technical-basics-and-calculation.html#c4622
% % Packing factor
% if NS<12
%     p = 1.25;
% elseif NS<16
%     p = 1.26;
% elseif NS<25
%     p = 1.27;
% elseif NS>=25
%     p = 1.28; % for 25<n<400
% end
% OD = p*sqrt(NS)*di_strand+ (increase in diameter by yarn serving)
% yarn = 0; % This is the layer that is around the whole cable. Here it is assumed that it is unserved, thus yarn =0
% di_strand = (di-yarn)/(p*sqrt(NS));
R_DC_single = CoilResistDC(l,di_strand);
% Rs_single = (4*1.7e-8/(di_strand^2))*(l/pi); % Total DC resistance of the cable
% Litz wire
R_DC = R_DC_single.*(1.015.^NB).*(1.025.^NC)./NS; % 
end