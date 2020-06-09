clear
clc
%% Read the information of belt data
% Import radar File
filename = '\\tsclient\d\Project\subject1_Normal_using_phone_2.07m_sitting_run2.csv';
[num,txt,raw]=xlsread(filename);
% [row, col]=size(num);
% radar_sample_num = (row+1)/2;

% Remove the NaN in data
num_time = num(:,1);
num_distance = num(:,2);
valid_radar_sample_time = num_time(~isnan(num_time));
valid_radar_sample_distance = num_distance(~isnan(num_distance));
valid_radar_sample = [valid_radar_sample_time';valid_radar_sample_distance']';
% Visualization
plot(valid_radar_sample_time,valid_radar_sample_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Applying SG Filter to remove noise from the data
rd = 2;  %Order of the SG Filter
fl = 101;  %Frame Length of the SG Filter
smtlb = sgolayfilt(valid_radar_sample_distance,rd,fl); %Applying the SG Filter to the Data file
subplot(2,1,1)
plot(valid_radar_sample_time, valid_radar_sample_distance);
title('Original')
grid

subplot(2,1,2)
plot(valid_radar_sample_time,smtlb);
title('Filtered')
grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time Frquency using Spectrogram which is STFT with default window size
Fs= 17;    %frequency of the Radar
figure;
spectrogram(valid_radar_sample_distance,[],[],[],Fs, 'yaxis');
%Time Frquency using Spectrogram which is STFT with reduced window
figure;
spectrogram(valid_radar_sample_distance,100,50,100,Fs, 'yaxis'); %after fine tuning window size=100,
%and noverlap=50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time Frequency using Continuous Wavelet Transformation
figure;
cwt(valid_radar_sample_distance,Fs);
%Time Frequency using Continuous Wavelet Transformation with FIne Scale
%Analysis
No= 7; %Number of Octaves
Nv=48; %Voices per Octaves
figure;
cwt(valid_radar_sample_distance,Fs,'NumOctaves',No,'VoicesPerOctave',Nv); %Applying Continuous Wavelet Transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparing CWT Results with the original signal and plotting
[cfs,frq]= cwt(valid_radar_sample_distance,Fs,'NumOctaves',No,'VoicesPerOctave',Nv);
tms= (0:numel(valid_radar_sample_distance)-1)/Fs; 
figure
subplot(2,1,1)
plot(tms,valid_radar_sample_distance)
axis tight
title('Signal and Scalogram')
xlabel('Time (s)')
ylabel('Distance')
subplot(2,1,2)
surface(tms,frq,abs(cfs))
axis tight
shading flat
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'yscale','log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



