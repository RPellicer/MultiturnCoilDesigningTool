function visuNoiseSweep(ProbeOpt)
    [a,b,c]=size(ProbeOpt);
        for k = 1:a
            for l = 1:b
                OptNoisesA(k,l,:) = abs(mean(ProbeOpt(k,l,1).Sens.SensPerNoise,2));
                if size(ProbeOpt,3)>1
                    OptNoisesB(k,l,:) = abs(mean(ProbeOpt(k,l,2).Sens.SensPerNoise,2));
                end
                OptSweepNr(l) = ProbeOpt(1,l,1).(ProbeOpt(1,l,1).Freqs.sweepPropA).(ProbeOpt(1,l,1).Freqs.sweepPropB);
            end
        end
    h = figure;
    ind = 1:size(OptNoisesA,2);
    loglog(OptSweepNr,squeeze(mean(OptNoisesA(:,:,:),1)),'-*')
    grid on
    legend(ProbeOpt(1,1,1).Sens.NameNoises)
    xlabel(ProbeOpt(1,1,1).Freqs.sweepPropB)
    ylabel('Sensitivity limit (T/sqrt(Hz))')
    title(['Sensitivity limit of each noise component sweep for different ' (ProbeOpt(1,1,1).Freqs.sweepPropB)])

    if size(ProbeOpt,3)>1
        hold on
        loglog(OptSweepNr,squeeze(mean(OptNoisesB(:,:,:),1)),'-.o')
        legend([ProbeOpt(1,1,1).Sens.NameNoises;...
            ProbeOpt(1,1,2).Sens.NameNoises])
        hold off
    end
end