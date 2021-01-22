%%
%     N_v = x(1); % Loops per layer
%     N_l = x(2); % Layers
%     do = x(3);  % [METERS] wire outer diameter insulation included
%     di = x(4);  % [METERS] cooper outer diameter
function stotal = calcSensTuned(amp,coil,matchNetw,freqs,Visual_glob)
    % Calculate
    stotal = sensTunedLinv(amp,coil,matchNetw,freqs,Visual_glob);    
end