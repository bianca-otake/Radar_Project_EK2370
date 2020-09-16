function [times, velocities] = CW_cont(data)
%CW_RADAR_ALGORITHM Summary of this function goes here
%   Detailed explanation goes here
fs = 44100;
Tp = 0.2;
zero_padding = 4;
ms_clutter = 0;
f=2.43*10^9;% 2.4 GHz


%constant
c = 3*10^8;   %m/s
N = Tp*fs;  %fs = sampling frequency and N = # of samples in time Tp
M_sweep = round(length(data(:,1))/N);

%Time domain matrix containing all the sampled value (velocities)
mat_time = reshape(data(1:(M_sweep-1)*(N),1),[N,M_sweep-1])';% inializing the matrix in time domain

%Frequency domain matrix
mat_freq = 20*log10(abs(fft(mat_time,zero_padding*N,2)));

matFby2 = mat_freq(:,1:length(mat_freq)/2);

max_value_freq = max(max(matFby2));
matFby2 = matFby2 - max_value_freq;


fd_max = fs/2;
fd_resolution = fd_max/length(matFby2);

v_max = c*fd_max/(2*f);
v_res = c*fd_resolution/(2*f);

%Axis
time_array = linspace(0,Tp*M_sweep,M_sweep);
v_array = linspace(0,v_max,length(matFby2));
% figure(2);
% imagesc(v_array,time_array,matFby2,[-45,0]);
% xlabel('velocity(m/s)');
% xlim([0 30])
% ylabel('Time(s)')
% colorbar
[~,I] = max(matFby2');
velocities = v_array(I);
velocities = velocities(velocities<200)
times = time_array(velocities<200)
end

