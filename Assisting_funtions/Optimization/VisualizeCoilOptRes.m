load('optimResults.mat')
disp(['Total computing TIME = ' num2str(floor(elapsedTime/(60*60))) ' hours, ' num2str(floor(elapsedTime/60)-60*floor(elapsedTime/(60*60))) ' mins, ' num2str(round(mod(elapsedTime,60))) ' secs' ]);
OptSens = zeros(size(ProbeOpt));
flagN = 0;
flagT = 0;
tempLitzStrNr = [];

[a,b,c]=size(ProbeOpt);
for k = 1:a
    for l = 1:b
        for m = 1:c
           OptSens(k,l,m) = ProbeOpt(k,l,m).Sens.SensScore;
           try
                OptSensPerNoise(k,l,m,:) = mean(ProbeOpt(k,l,m).Sens.SensPerNoise,2);
           catch
           end
           OptType(k,l,m) = ProbeOpt(k,l,m).Type(1);
           OptSweepNr(k,l,m) = ProbeOpt(k,l,m).(ProbeOpt(k,l,m).Freqs.sweepPropA).(ProbeOpt(k,l,m).Freqs.sweepPropB); %LitzStrandNr_glob;
        end
    end
end
%% Generate the list
%  ProbeOpt(count,litzCount,countType)
OptSensT2 = reshape(OptSens,a,b*c,1);
OptTypeT2 = squeeze(reshape(OptType(1,:,:),1,b*c,1));
OptSweepNrT2 = squeeze(reshape(OptSweepNr(1,:,:),1,b*c,1));
for i = 1: (b*c)
    list(i) = {[OptTypeT2(i) num2str(abs(OptSweepNrT2(i)))]};
end
% list(regexp(list,'.'))=[]
list = regexprep(list,'[.]','');
%% do the statistical analysis
maskS = ~~OptSensT2;
meanOptSensT = sum(OptSensT2,1)./sum(maskS,1);
%     meanOptSensT = mean(OptSensT,1);
stdOptSensT = sqrt(1./sum(maskS,1).*sum(maskS.*(OptSensT2-repmat(meanOptSensT,size(OptSensT2,1),1)).^2));
%     stdOptSensT = std(OptSensT,[],1);
minOptSensT = min(OptSensT2+~maskS*1000,[],1); % Adds a big value to zero values so they don't show as minimum
try
    array2table(OptSensT2, 'VariableNames', list)
catch
end
figure(16)
semilogy(OptSensT2,'-*')
grid on
legend(list), xlabel('Iteration Nr.'), ylabel('Sensitivity (T/sqrt(Hz))')
figure(17)
% bar(meanOptSensT2,'grouped')
bar(meanOptSensT)
colormap(summer)
hold on
plot(minOptSensT,'*r')
set(gca, 'XTick',1:length(list));
set(gca, 'XTickLabel',list),
set(gca, 'XTickLabelRotation',45);
% legend(listNums)
ylabel('Sensitivity (T/sqrt(Hz))')
title('Sensitivity with average, std, and best results')
errorb(meanOptSensT(1,:),stdOptSensT(1,:));

% Compare results in parallel
if size(ProbeOpt,2)>1
    figure(18)
    loglog(OptSweepNr(1,:,1),(squeeze(mean(OptSens,1))),'-*')
%     set(gca, 'XTick',1:length(OptSweepNr(1,:,1)));
%     set(gca, 'XTickLabel',OptSweepNr(1,:,1)),
    set(gca, 'XTickLabelRotation',45);
    xlabel(ProbeOpt(1,1,1).Freqs.sweepPropB)
    ylabel('Sensitivity (T/sqrt(Hz))')
    title(['Sensitivity sweep for different ' (ProbeOpt(1,1,1).Freqs.sweepPropB)])
    if size(ProbeOpt,3)>1
        legend(ProbeOpt(1,1,1).Type,ProbeOpt(1,1,2).Type)
    else
        legend(ProbeOpt(1,1,1).Type)
    end
    grid on
    
    figure(19)
%     values = squeeze(mean(OptSens,1));
    values = squeeze(min(OptSens,[],1));
    maxValuesT = max(values);
    maxValues = ones(1,size(values,1))'*maxValuesT;
    normValues = values./maxValues;
    loglog(OptSweepNr(1,:,1),normValues,'-*')
    set(gca, 'XTickLabelRotation',45);
    xlabel(ProbeOpt(1,1,1).Freqs.sweepPropB)
    ylabel('Sensitivity difference')
    title(['Sensitivity difference sweep for different ' (ProbeOpt(1,1,1).Freqs.sweepPropB)])
    if size(ProbeOpt,3)>1
        legend(ProbeOpt(1,1,1).Type,ProbeOpt(1,1,2).Type)
    else
        legend(ProbeOpt(1,1,1).Type)
    end
    grid on

end
%% Visualize the noise for the swept variable
if maxSweepCount>1
    visuNoiseSweep(ProbeOpt)
end
