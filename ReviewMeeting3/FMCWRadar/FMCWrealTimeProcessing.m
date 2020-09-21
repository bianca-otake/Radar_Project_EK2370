%% Real  time CW audio processing
clc, clear all, close all

% initialize data arrays
recordedData = [0 0]
t_array = [0]
r_array = [0]

% create recording object
Fs = 44100 ; 
nBits = 16 ; 
nChannels = 2 ; 
ID = -1; % default audio input device 
recObj = audiorecorder(Fs,nBits,nChannels,ID);

Tp = 1 % seconds
f = 0
h = figure(1);
t = 0
while f == 0 
    t = t + Tp 
    recordblocking(recObj,Tp);
    audioData = recObj.getaudiodata;
    recordedData = [recordedData;audioData];

    [time, ranges] = FMCW_range(recordedData);
%     t_array = [t_array;time' + t_array(end)]
%     r_array = [r_array; ranges'];
%     recordedData = [recordedData;audioData];
    figure(1)
%     plot(t_array,r_array)
    plot(time,ranges)

    isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
    if isKeyPressed
     break
    end
end

figure(2)
plot(recordedData(:,1))
title('Sampled Signal')
xlabel('Sample number [-]')

%% plot the total data and its spectrogram
FMCW_spectrogram_range(recordedData)
% [time, ranges] = FMCW_range(recordedData);
figure(4);
% plot(t_array,r_array)


