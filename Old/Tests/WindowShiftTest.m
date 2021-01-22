close all, clear all, clc
fs = 1e4; % 1kHz
dt = 1/fs; 
windowLength = 0.1; % seconds
t = (0:dt:windowLength-dt);
res = 1/windowLength;
disp(['Resolution limit = ',num2str(res),' Hz'])
nrf_fft = 20000;


% Frequencies that will compose the signal (in Hertzs)
f = [80, 100, 102, 110, 120, 140]; % Radians

diff_f = f(2:end)-f(1:end-1);
reqRes = min(diff_f);
disp(['Min. required resolution = ',num2str(reqRes),' Hz'])

s_indiv = sin(f'*(2*pi)*t);
s = sum (s_indiv,1);

% Visualize the challenge
figure(1)
plot(t,s_indiv);
hold on
plot(t,s,'k','LineWidth',3)

%% FFT
Yf = (fft(s,nrf_fft,2)/length(s))*windowLength;
fAxistemp = fs*linspace(0,1,nrf_fft+1);
fAxis = fAxistemp(1:end-1);
% Visualize
figure(2);
CH = abs(Yf);
plot(fAxis,CH,'r');

%% Manual rearrange of windows
acquiredPeriodNr = windowLength*f;
optimumPeriodNr = floor(acquiredPeriodNr);
shift = ((acquiredPeriodNr - optimumPeriodNr)./f)/dt;

%% Add windowing sequence
concNr = 1/(reqRes*windowLength);
% concNr = 100;
s_long = repmat(s,1,concNr);
t_long = t;
for i=2:concNr
    t_long = [t_long (t+t_long(end)+dt)];
end;
figure(1),
plot(t_long, s_long,'-.b','LineWidth',1.5)

%% FFT
Yf2 = (fft(s_long,nrf_fft,2)/length(s_long))*windowLength*concNr;
fAxistemp2 = fs*linspace(0,1,nrf_fft+1);
fAxis2 = fAxistemp2(1:end-1);
% Visualize
figure(3);
CH2 = abs(Yf2);
plot(fAxis2,CH2,'r');

% Experimenting
%% Add windowing sequence
concNr = 1/(reqRes*windowLength);
% concNr = 100;
currentshift = round(shift(3));
s_temp = s(1:end-currentshift);
s_long = repmat(s_temp,1,concNr);
t_temp = t(1:end-currentshift);
t_long = t_temp;
for i=2:concNr
    t_long = [t_long (t_temp+t_long(end)+dt)];
end;
%% FFT
Yf4 = (fft(s_long,nrf_fft,2)/length(s_long))*windowLength*concNr;
fAxistemp4 = fs*linspace(0,1,nrf_fft+1);
fAxis4 = fAxistemp4(1:end-1);
% Visualize
figure(6);
CH4 = abs(Yf4);
plot(fAxis4,CH4,'r');

% Visualyze the problem 
s3_long = repmat(s_indiv(3,:),1,concNr);
figure(4),plot(s3_long);
Yf3 = (fft(s3_long,nrf_fft,2)/length(s3_long))*windowLength*concNr;
% Visualize
figure(5);
CH3 = abs(Yf3);
plot(fAxis2,CH3,'r');