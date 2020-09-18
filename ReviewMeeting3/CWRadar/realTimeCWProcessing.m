%% Real  time CW audio processing
clc, clear all, close all

recordedData = []
t_array = [0]
v_array = [0]

% Set for how long the live processing should last
endTime = 20;
audioFrameLength = 3200;
Fs = 44100 ; 
nBits = 16 ; 
nChannels = 2 ; 
ID = -1; % default audio input device 
Tp = 0.3 % 0.3 seconds 
recObj = audiorecorder(Fs,nBits,nChannels,ID);
f = 0
h = figure(1);
t = 0
while f == 0 
    t = t + Tp 
    recordblocking(recObj,1);
    audioData = recObj.getaudiodata;
    
    [time, velocities ] = CW_cont(audioData);
    t_array = [t_array;time' + t_array(end)]
    v_array = [v_array; velocities'];
    recordedData = [recordedData;audioData(:,1)];
%     figure(2);
%     plot(audioData(:,1))
    figure(1)
    plot(t_array,v_array)
    isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
    if isKeyPressed
     break
    end
end
figure(2)
plot(recordedData)
title('Sampled Signal')
xlabel('Sample number [-]')

CW_spectrogram(recordedData)


