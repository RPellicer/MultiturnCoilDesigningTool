function Probe = getProbeDetails(tProbe,graphs)
    try
        disp(['Optimization TIME = ' num2str(floor(tProbe.Sens.OptProcTime/(60*60))) ' hours, ' num2str(floor(tProbe.Sens.OptProcTime/60)-60*floor(tProbe.Sens.OptProcTime/(60*60))) ' mins, ' num2str(round(mod(tProbe.Sens.OptProcTime,60))) ' secs' ]);
    end
    if ~isfield(tProbe.Freqs,'f0_Offset')
        tProbe.Freqs.f0_Offset = 0;
    end
    disp('#######################################################################################################################################################')
    switch tProbe.Type 
        case 'Tuned'
            if tProbe.Coil.Coil_heigth_cent == 1
                L = abs(tProbe.Coil.hBase_ROI-tProbe.Coil.hfinal_ROI);
                shift = (tProbe.Coil.do*tProbe.Coil.N_v)/2;
                tProbe.Coil.hBase_ROI = -L/2+shift;         % Closes part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
                tProbe.Coil.hfinal_ROI = L/2+shift;%tProbe.Coil.Rout_glob; % Furthest part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
            end
            [~, Probe] = sensTunedLinv(tProbe,graphs);
            if tProbe.Coil.Rin ~= tProbe.Coil.Rin_glob
                warning('5 getProbeDetails: Rin different to Rin differes from Rin_glob') 
            end
            disp(['Tuned probe, StrandNr = ' num2str(Probe.Coil.LitzStrandNr_glob)])
            disp(['Estimated optimal coil has ' num2str(Probe.Coil.N_v) ' loops per layer, ' num2str(Probe.Coil.N_l) ' layers, with spacing of  ' num2str(Probe.Coil.do) ' in between wires and cooper diameter of ' num2str(Probe.Coil.di) ' (all in meters)']);
            disp(['C1 = ' num2str(Probe.MatchNetw.C1) ', C2 = ' num2str(Probe.MatchNetw.C2)])
            disp(['Peak sensitivity = ' num2str(Probe.Sens.Sens_peak)  ', Average sensitivity = ' num2str(Probe.Sens.SensMean)])
            disp('-------------------------------------------------------------------------------------------------------------------------------------------------------')
            % Data to use with Pspice simulator to compare results:
            disp('Data to use into Pspice simulator to double check results:')
            disp(['Coil: Rs = ' num2str(Probe.Coil.Rs0) ', Ls = ' num2str(Probe.Coil.Ls) ', Cs = ' num2str(Probe.Coil.Cs) ', PspiceFreqDep = ' num2str(Probe.Coil.PspiceFreqDep)])
            disp(['Matching network: C1 = ' num2str(Probe.MatchNetw.C1) ', C2 = ' num2str(Probe.MatchNetw.C2), ', C1_ESR = ' num2str(Probe.MatchNetw.C1_ESR) ', C2_ESR = ' num2str(Probe.MatchNetw.C2_ESR)])
            disp(['Rmatch = ' num2str(Probe.MatchNetw.Rmatch) ', and self-resonant freq = ' num2str(Probe.Coil.SelfRes)])
        case 'Non-tuned'
            if tProbe.Coil.Coil_heigth_cent == 1
                L = abs(tProbe.Coil.hBase_ROI-tProbe.Coil.hfinal_ROI);
                shift = (tProbe.Coil.do*tProbe.Coil.N_v)/2;
                tProbe.Coil.hBase_ROI = -L/2+shift;         % Closes part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
                tProbe.Coil.hfinal_ROI = L/2+shift;%tProbe.Coil.Rout_glob; % Furthest part of the ROI/ Not needed if tProbe.Coil.emfNear = 0;
            end
            [~, Probe] = sensitivityI2V(tProbe,graphs);
            disp(['Non-tuned probe, StrandNr = ' num2str(Probe.Coil.LitzStrandNr_glob)]);
            disp(['Estimated optimal coil has ' num2str(Probe.Coil.N_v) ' loops per layer, ' num2str(Probe.Coil.N_l) ' layers, with spacing of  ' num2str(Probe.Coil.do) ' in between wires and cooper diameter of ' num2str(Probe.Coil.di) ' (all in meters)']);
            disp(['Peak sensitivity = ' num2str(Probe.Sens.Sens_peak)  ', Average sensitivity = ' num2str(Probe.Sens.SensMean)])
            % Data to use with Pspice simulator to compare results:
            disp('-------------------------------------------------------------------------------------------------------------------------------------------------------')
            disp('Data to use into Pspice simulator to double check results:')
            disp(['Coil: Rs = ' num2str(Probe.Coil.Rs0) ', Ls = ' num2str(Probe.Coil.Ls) ', Cs = ' num2str(Probe.Coil.Cs) ', PspiceFreqDep = ' num2str(Probe.Coil.PspiceFreqDep)])     
            disp(['Self-resonant freq = ' num2str(Probe.Coil.SelfRes)])
        otherwise
            warning('Unexpected coil type. It is not "Tuned" nor "Non-tuned".')
    end
    if Probe.Sens.PenalyzeScaling ~= 1
        disp(['The optimization was aproaching the resonant frequency Fres = ' num2str(Probe.Coil.SelfRes ) ', f0 = '  num2str(Probe.Freqs.f0_glob)])
    end
    disp('#######################################################################################################################################################')
    if Probe.Coil.PenalyzeScaling > 1
        disp(['It has been penalized by ' num2str(Probe.Coil.PenalyzeScaling)])
    end
    disp(' ')
end