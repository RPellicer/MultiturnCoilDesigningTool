%% This code reads the raw signal recorded from Labview of the sensor and 
% the Bp current, syncronizes the signals looking at the current and
% averages them.
% Centre for Advanced Imaging, UQ, Brisbane, Ruben Pellicer 5/10/2015

clear all, close all, clc
%% Extracting the data from the file and the structure
% Select the file GUI
[filename, PathName, dataStruct] = selectLabviewFile();
% Get data from the file
data = getLabviewdata(dataStruct);
% Syncronize the signals with the current of Bp
ringTime = 0.02; % seconds from the plot
data.fs = 64000; % Overwrite because it is wrong in the file
[data.signalSync, data.sampleSyncLength] = syncBpcurrent(data, ringTime);
% Remove offset (DC)
data.signalAC = data.signalSync - repmat(mean(data.signalSync,2),1,size(data.signalSync,2));
% FFT
[data.v_f, data.freq] = v_fft(data.signalAC,data.fs,2); % FOLDED Spectrum
% Averaged FFT
timeAvedata = mean(data.signalAC,1);
% figure = h;
% plot(data.time(1:data.sampleSyncLength),timeAvedata)
[data.v_f, data.freq] = v_fft(timeAvedata,data.fs,2); % FOLDED Spectrum
% RMS in frequency domain
V_RMS_num = RMSspectrum(data.v_f,data.freq);
% RMS in time domain
V_RMS = rms(data.signalAC,2);

%% Visualize spectrum
% Rectangular window FFT style
tracks = [29 3 50 3];
pic = figure;
signalNr = 1;
loglog(data.freq,abs(data.v_f(signalNr,:)),'b')
loglog(data.freq,abs(data.v_f(signalNr,:)),'b')
hold on, grid on
% Visualizing with pwelch for smoothness
sig = data.signalAC(signalNr,:);
% Hanning window without overlap
nx = max(size(sig));
na = 16; % Number of averages
w = hanning(floor(nx/na));
[Pxx,f]=pwelch(sig,w,0,[],data.fs);
loglog(f(tracks(2):tracks(3)),sqrt(Pxx(tracks(2):tracks(3))),'m')
% Hanning window with overlap
na = 100; % Number of averages
w = hanning(floor(nx/na));
[Pxx2,f2]=pwelch(sig,w,[],[],data.fs);
loglog(f2(tracks(4):end),sqrt(Pxx2(tracks(4):end)),'r')
xlabel('freq. (Hz)'), ylabel('Output voltage noise V/Hz^1^/^2')
title('Output voltage noise spectral density')
legend('Measured low freq.','Measured medium freq.','Measured high freq.')
