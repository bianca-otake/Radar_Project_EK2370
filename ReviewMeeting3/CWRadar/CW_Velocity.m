close all
clear all
clc
 
load('velocity_data.mat') % to read file
velocity=-1*data(:,1); % Inversion due to Amplifier gain is negative
N_total=length(velocity); % Total number of samples in read file
T=N_total/fs; % Total time duration of signals in read file
fs = 44100;
Tp=0.3; % Tp in seconds
% fs is the sampling frequency in HZ
N= Tp*fs; % Number of sample within Tp.
M_sweeps=round(T/Tp); % Number of rows in time domain matrix

mat_time=reshape(velocity(1:(M_sweeps-1)*(N)),[M_sweeps-1,N]);% inializing the matrix in time domain
% for i = 1:M_sweeps-1
%     mat_time(i,:) = velocity(1 + N*(i-1):i*N); % assiging amplitude
% end
%% Performing the MS(Mean Substract) clutter rejection
Ms_mean=mean(velocity); % Mean of whole matrix
mat_time=mat_time-Ms_mean; % MS clutter removal
%% Performing zero padding in time domain matrix
mat_time_padding = [zeros(M_sweeps,N*4) mat_time zeros(M_sweeps,N*4)];
% for i = 1:M_sweeps-1
%     mat_freq(i,:)= 10*log(abs((fft(mat_time(i,:),4*N))));
%     %mat_freq(i,:)= (((fft(mat_time(i,:),4*N))));
%  end

% mat_freq= (abs((fft(mat_time,4*N))));
%mat_freq= ((fft(mat_time_padding)));
mat_freq = 10*log(abs(fft(mat_time_padding)));
max_value = max(max(mat_freq)); %max value of the matrix
mat_freq =mat_freq-max_value; %normaliz
time_array(:,1) = 0:Tp:Tp*M_sweeps;
F = fs/2;
N_freq=length(mat_freq)/2;
%velocity_mat = zeros(M_sweeps,N_freq/2);
for i = 1:M_sweeps
    mat_freq_modify(i,:) = mat_freq(1 + N_freq*(i-1):i*N_freq);
end
f=2.43*10^9;% 2.4 GHz

for i=1:M_sweeps
    dopper(i)=max(mat_freq_modify(i,:));
    velocity_measured(i,1)=(3*10^8/(2*f))*dopper(i);
end

%imagesc(velocity_measured,time_array,abs(fft(mat_time_padding)))
imagesc(velocity_measured,time_array,abs(fft(mat_time)))
xlabel('velocity(m/s)')
ylabel('Time(s)')
colorbar