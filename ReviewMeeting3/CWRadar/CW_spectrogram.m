function [times, velocities] = CW_spectrogram(data)
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
figure(100);
imagesc(v_array,time_array,matFby2,[-45,0]);
xlabel('velocity(m/s)');
xlim([0 30])
ylabel('Time(s)')
colorbar
%% 

    Nmax =3
    [Tp,~] = islocalmax(matFby2,2,'MaxNumExtrema',Nmax);
    ranges = zeros(size(Tp,1),Nmax);
    for i = 1:size(Tp,1)
       ranges(i,:) = range_array(Tp(i,:)) ;
    end
    figure()
    subplot(2,1,1);
    hold on
    for i = 1:Nmax 
        scatter(time_array,ranges(:,i),'x')
    end
    xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('FMCW Radar Local and global maximum','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    ylim_inf = 0;
    ylim_sup = 100;
    ylim([ylim_inf ylim_sup])  
    % global max
    [~,I] = max(matFby2');
    velocities = v_array(I);
    velocities = velocities(velocities<200)
    times = time_array(velocities<200)
    
    subplot(2,1,2);
    plot(times,velocities)
    xlabel('Velocity(m/s)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    set(gca,'FontSize',10,'FontWeight','bold');
end

