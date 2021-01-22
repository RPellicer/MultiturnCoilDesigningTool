% AmplMode=> 
%  0 = Non-tuned
%  1 = Tuned
function [s,Probe] = CalcSensFromDesign(NrLoopsPerLayer,NrLayers,OutWireDiam,CoopDiam, AmplMode)
%     run'C:\Users\s4370509\Documents\MATLAB\MultiturnDetector\Optimization\setGlobalVariables.m')
%     run('setGlobalVariables_MiguelField.m')
%     run('setGlobalVariables_MiguelNMR.m')
    run('setGlobalVariables.m')
    visu =1;
    tProbe.Coil.N_v = NrLoopsPerLayer;
    tProbe.Coil.N_l = NrLayers;
    tProbe.Coil.do = OutWireDiam;    
    tProbe.Coil.di = CoopDiam;
    % Set other needed variables
    tProbe.Coil.Rin = tProbe.Coil.Rout_glob - tProbe.Coil.N_l*tProbe.Coil.do;
    if tProbe.Coil.emfNear == 1
        tProbe.Coil.Rin = tProbe.Coil.Rin_glob;
    end
    percInsul = 0.99;
    %     Centre the coil at origing if flag "tProbe.Coil.Coil_heigth_cent" is on 
    if tProbe.Coil.Coil_heigth_cent == 1
        L = abs(tProbe.Coil.hBase_ROI-tProbe.Coil.hfinal_ROI);
        shift = (tProbe.Coil.do*tProbe.Coil.N_v)/2;
        tProbe.Coil.hBase_ROI = -L/2+shift;         % Closes part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
        tProbe.Coil.hfinal_ROI = L/2+shift;%tProbe.Coil.Rout_glob; % Furthest part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
    end
    switch AmplMode
        case 0 % Non-tuned
            if (tProbe.Coil.do*percInsul) <= tProbe.Coil.di
                tProbe.Coil.PenalyzeScaling = tProbe.Coil.PenalyzeScaling*(1+ 5*((1/(1-percInsul))*((tProbe.Coil.di-tProbe.Coil.do*percInsul)/tProbe.Coil.di))); % The second part makes the sensitivity 5 times worse when do=di. It can be adjusted
            end
            [s,Probe] = sensitivityI2V(tProbe,visu);
        otherwise % Tuned
            if (tProbe.Coil.do*percInsul) <= tProbe.Coil.di
                tProbe.Coil.PenalyzeScaling = tProbe.Coil.PenalyzeScaling*(1+ 5*((1/(1-percInsul))*((tProbe.Coil.di-tProbe.Coil.do*percInsul)/tProbe.Coil.di))); % The second part makes the sensitivity 5 times worse when do=di. It can be adjusted
            end
            [s,Probe] = sensTunedLinv(tProbe,visu);
    end
end