function [sensScore, probe] = sensTunedLinvManual(N_v, N_l, do, di)
    close all
    setGlobalVariables()
    tProbe.Coil.N_v = N_v; % Loops per layer
    tProbe.Coil.N_l = N_l; % Layers
    tProbe.Coil.do = do;  % [METERS] wire outer diameter insulation included
    tProbe.Coil.di = di;  % [METERS] cooper outer diameter
    tProbe.Coil.Rin = tProbe.Coil.Rout_glob - tProbe.Coil.N_l*tProbe.Coil.do; % Coil internal radius    f0 = (max(freqs_glob)+min(freqs_glob))/2;
    visu = 1;
    [sensScore, probe] = sensTunedLinv(tProbe,visu);
end