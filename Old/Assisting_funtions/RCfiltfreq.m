%% Calculcating resonant frequency of a filter
% Enter the known values of two of R, C or the Frequency and the function
% will find the missing one. e.g.:
% [R,~,~] = RCfiltfreq([],80e-9,200)
function  [R,C,F] = RCfiltfreq(R,C,F)
    if isempty(R)
      R = 1/(2*pi*F*C);
    elseif isempty(C)
      C = 1/(2*pi*R*F);
    elseif isempty(F)
      F = 1/(2*pi*R*C);
    end
end