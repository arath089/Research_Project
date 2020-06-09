%% Cheynes-Stokes Simulation
%% Parameters
fs=256; % Sampling Rate
br=15; % Breathing Rate
cycle=50; % Cheyne-Stokes cyle time in secods
apnea=10; % Apnea length
time=10; % Total time in minutes
VarEnv=0.1; % Variance of envelope noise
VarPhase=0.1; % Variance of phase noise in breathing signal

%% Generate Signal
t=linspace(0,time*60, time*60*fs); % Time vector
Sn=cos(2*pi*(br/60).*t+VarPhase*randn(1,time*60*fs)); % Normal Breathing signal
% Construct hyperpnea/apnea window function 
win=blackman(cycle*fs)'; % crescendo-diminuendo pattern represented by Blackman window
apneawin=zeros(1,apnea*fs); % apnea window represented by zeros
envelope=[win,apneawin]; % One period of envelope cycle
% Periodically extend window function
while length(envelope)<length(Sn)
    envelope=[envelope,envelope];
end
envelope=envelope(1:length(Sn))+VarEnv.*randn(1,time*60*fs); % Trim length of window to match signal length
CSB=Sn.*envelope; % Cheynes Stokes breathing pattern


%% Plot Signals
subplot(2,1,1)
plot(t/60,Sn)
title('Normal Breathing')
xlabel('time (min)')
ylabel ('amplitude')
subplot(2,1,2)
plot(t/60,CSB)
title('Cheyne-Stokes Breathing')
xlabel('time (min)')
ylabel ('amplitude')
