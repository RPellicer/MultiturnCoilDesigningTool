%% Quick coil modulation principle test
fs = 1000;
dt=1/fs;
acqWindow = 1;
t = 0:dt:acqWindow-dt;
signalAmplitude = 1;
fsample = 5;
fcarr = 100;
carrmax = 1; %0.9;
carrmin = 0; %0.3;
feqRes = 1/acqWindow;
signalSample = signalAmplitude*cos(2*pi*fsample*t);
signalCarrier = (carrmax+carrmin)/2+((carrmax-carrmin)/2)*cos(2*pi*fcarr*t);
signalCoil2 = signalSample.*signalCarrier;

subplot(5,1,1);
plot(t,signalSample,'b'),
title('Signal coming from the sample'),
xlabel('Time (s)'),ylabel('Amplitude(u)'),
subplot(5,1,2);
plot(t,signalCarrier,'k'),
title('Carrier signal'),
xlabel('Time (s)'),ylabel('Amplitude(u)'),
subplot(5,1,3);
plot(t,signalCoil2,'r'),
title('Modulated signal'),
xlabel('Time (s)'),ylabel('Amplitude(u)'),

%% FFT
Y_signalSample = fft(signalSample,[],2)/acqWindow;
Y_coilpos = fft(signalCarrier,[],2)/acqWindow;
Y_fluxCoil = fft(signalCoil2,[],2)/acqWindow;
Yf_signalSample = Y_signalSample;
Yf_coilpos = Y_coilpos;
Yf_fluxCoil = Y_fluxCoil;
% //Calculate frequency axis
NFFT = size(Y_signalSample,2);
df = fs/NFFT;
fAxisTemp = fs*linspace(0,1,NFFT+1);
fAxis = fAxisTemp(1:end-1);
% Visualize
subplot(5,1,4);
CH_signalSample = abs(Yf_signalSample);
CH_coilpos = abs(Yf_coilpos);
CH_fluxCoil = abs(Yf_fluxCoil);
plot(fAxis,CH_signalSample,'b'), hold on
plot(fAxis,CH_coilpos,'k'),
plot(fAxis,CH_fluxCoil,'r'),
hold off
title('Spectra of signal coming from the sample, coil rotating frequency and resultant flux accross the coil'),
xlabel('Frequency (Hz)'),ylabel('Amplitude(u)'),
legend(['Spin precesion (Hz) = ',num2str(fsample)],['Carrier (Hz) = ',num2str(fcarr)],['Flux experienced at the coil (Hz)']);

%% Enhancement
sensitivity = fAxis/fsample;
sensCH_signalSample = CH_signalSample.*sensitivity;
sensCH_coilpos = CH_coilpos.*sensitivity;
sensCH_fluxCoil = CH_fluxCoil.*sensitivity;
% figure(2);
subplot(5,1,5);
half = floor(length(fAxis)/2);
plot(fAxis(1:half),sensCH_signalSample(1:half),'b'), hold on
plot(fAxis(1:half),sensCH_coilpos(1:half),'k'),
plot(fAxis(1:half),sensCH_fluxCoil(1:half),'r'),
hold off
title('Spectra of Voltage mesured by the coil (v=dflux/dt taken into account)'),
xlabel('Frequency (Hz)'),ylabel('Amplitude(u)'),
legend(['Spin precesion (Hz) = ',num2str(fsample)],['Carrier (Hz) = ',num2str(fcarr)],['Flux experienced at the coil (Hz)']);

%% Estimation
Bp = 0.1;
Bm = 1082*54e-6;
Carrier = 5e9;
B0comp = 1.5;
fBm = Bm * 127e6/3;
fcomp = B0comp * 127e6/3;
eff = 0.5*0.5*0.4;
SNRcomp = B0comp * fcomp;
SNRulf = Bp * (Carrier+fBm) * eff;
disp(SNRulf/SNRcomp)


Bp = 1.5;
Bm = 1.5;
Carrier = 200e9;
B0comp = 1.5;
fBm = Bm * 127e6/3;
fcomp = B0comp * 127e6/3;
eff = 0.5*0.5*0.4;
SNRcomp = B0comp * fcomp;
SNRulf = Bp * (Carrier+fBm) * eff;
disp(SNRulf/SNRcomp)