function [signalSync, sampleSyncLength] = syncBpcurrent(data,ringTime)
    [~,Ind] = min(data.signal(:,:,2),[],2);
%     ringSample = ringTime*data.fs;
    ringSample = ringTime*data.fs; % This has to be removed
%     It uses the time after ringing minus the difference in trigering
    sampleSyncLength = data.sampleLength -ringSample-(max(Ind)-min(Ind));
    signalSync = zeros(data.sampleNr,sampleSyncLength);
    for i = 1: data.sampleNr
%         It takes the sample from the ringing time plus the difference
%         from its trigger and the maximum triger for the samplesynclength
        signalSync(i,:) = data.signal(i,[(ringSample+(Ind(i)-min(Ind))):((ringSample+(Ind(i)-min(Ind)))+sampleSyncLength-1)],1);
        b(i)=(ringSample+(Ind(i)-min(Ind)));
    end
    h = figure;
    plot(data.time(1:sampleSyncLength),squeeze(signalSync(:,:))),xlabel('Time(s)'), ylabel('Volts'),grid on % Visualize all the signals
end