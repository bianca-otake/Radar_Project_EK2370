function CW_radar_algorithm(fc_trial, Tp, zero_padding, ms_clutter)
%CW_RADAR_ALGORITHM Summary of this function goes here
%   Detailed explanation goes here
close all
load('velocity_data.mat') % to read file

figure1 = figure; hold on;
plot(-data(:,1),'b')
plot(data(:,2),'r')
xlabel('Data sample number','FontName','Times')
xlim([0.5*10^5 2.5*10^5])
ylabel('Amplitude','FontName','Times')
set(gca,'FontSize',10,'FontWeight','bold')
hold off

%constant
c = 3*10^8;   %m/s
N = floor(Tp*fs);  %fs = sampling frequency and N = # of samples in time Tp
M_sweep = round(length(data(:,1))/N);

%Time domain matrix containing all the sampled value (velocities)
mat_time = reshape(data(1:(M_sweep-1)*(N),1),[N,M_sweep-1])';% inializing the matrix in time domain

%MS clutter rejection
switch ms_clutter 
    case 1
        mean_time = mean2(mat_time);
        mat_time = mat_time - mean_time;
end

%Frequency domain matrix
mat_freq = 20*log10(abs(fft(mat_time,zero_padding*N,2)));

figure2 = figure;
plot(mat_freq(1,:))

matFby2 = mat_freq(:,1:length(mat_freq)/2);
figure3 = figure;
plot(matFby2(1,:))

max_value_freq = max(max(matFby2));
matFby2 = matFby2 - max_value_freq;
figure4 = figure;
plot(matFby2(1,:))

fd_max = fs/2;
fd_resolution = fd_max/length(matFby2);

v_max = c*fd_max/(2*fc_trial);
v_res = c*fd_resolution/(2*fc_trial);

%Axis
time_array = linspace(0,Tp*M_sweep,M_sweep);
v_array = linspace(0,v_max,length(matFby2));

figure5 = figure;
imagesc(v_array,time_array,matFby2,[-45,0]);
xlabel('{\it{velocity}}(m/s)','FontName','Times');
xlim([0 30])
ylabel('{\it{time}}(m/s)','FontName','Times')
title('Pulse time TP = ' + string(Tp) + 's, Center Freqency fc =' + string(fc_trial/(10^9)) + 'GHz','FontName','Times')
set(gca,'FontSize',10,'FontWeight','bold')
colorbar

%saveas(figure1,'OUTPUT\Signal.png')
%saveas(figure2,'OUTPUT\TP_variation\EntireRawDopplerFreq_0.2.png')
saveas(figure5,'OUTPUT\MS_clutterRejection\Output_' + string(fc_trial/(10^9)) + '.png')

end

