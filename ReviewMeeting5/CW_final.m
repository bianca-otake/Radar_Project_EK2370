close all
clear all
clc

Tp = 0.2;
fc_trial = 5.75*10^9;
ms_clutter = 1;
zero_padding = 4;
figure_no=3;

%CW_RADAR_ALGORITHM Summary of this function goes here

data = audioread('IQdemodulation.wav');
% data = audioread('multiple_target_CW_rev3.wav');

fs=100*10^3; % to read file
I = data(:,1);
Q = data(:,2);

data_2(:,1)= I+i*Q;
data_2(:,2)= I-i*Q;

% figure1 = figure; hold on;
% plot(-data(:,1),'b')
% %plot(data(:,2),'r')
% xlabel('Data sample number','FontName','Times')
% xlim([0.5*10^5 2.5*10^5])
% ylabel('Amplitude','FontName','Times')
% set(gca,'FontSize',10,'FontWeight','bold')
% hold off
for d = 1:2
    %constant
    c = 3*10^8;   %m/s
    N = floor(Tp*fs);  %fs = sampling frequency and N = # of samples in time Tp
    M_sweep = round(length(data_2(:,1))/N)-1;
    
    %Time domain matrix containing all the sampled value (velocities)
    mat_time = reshape(data_2(1:(M_sweep)*(N),d),[N,M_sweep])';% inializing the matrix in time domain
    
    %MS clutter rejection
    switch ms_clutter
        case 1
            mean_time = mean2(mat_time);
            mat_time = mat_time - mean_time;
    end
    
    %Frequency domain matrix
    mat_freq = 20*log10(abs(fft(mat_time,zero_padding*N,2)));
    
    % figure2 = figure;
    % plot(mat_freq(1,:))
    [m,n]=size(mat_freq);
    matFby2 = mat_freq(:,1:n/2);
    % figure3 = figure;
    % plot(matFby2(1,:))
    
    max_value_freq = max(max(matFby2));
    matFby2 = matFby2 - max_value_freq;
    % figure4 = figure;
    % plot(matFby2(1,:))
    
    fd_max = fs/2;
    fd_resolution = fd_max/length(matFby2);
    
    v_max = c*fd_max/(2*fc_trial);
    v_res = c*fd_resolution/(2*fc_trial);
    
    %Axis
    time_array = linspace(0,Tp*M_sweep,M_sweep);
    v_array = linspace(0,v_max,length(matFby2));
    
    % Axis for velocity vs time plot
    % [~,I] = max(matFby2');
    % velocities = v_array(I);
    
    A=matFby2';
    [aa,indices]=sort(A,'descend');
    velocities_2 = v_array(indices(1,:));
    %velocities_3 = v_array(indices(2,:));
    
        if (d == 1)
        print_f=sprintf('Moving Away');
        else
        print_f=sprintf('Moving Towards');
        end
    figure(3)
    subplot(1,2,d)
    imagesc(v_array,time_array,matFby2,[-45,0]);
    xlabel('Velocity(m/s)','FontName','Times');
    xlim([0 30])
    ylabel('Time(s)','FontName','Times')
%     if (d==1)
%         title('Moving away','Pulse time TP = ' + string(Tp) + 's, Center Freqency fc =' + string(fc_trial/(10^9)) + 'GHz, Zero Padding = '+ string(zero_padding) + 'N','FontName','Times')
%     else
%         title('Moving towards','Pulse time TP = ' + string(Tp) + 's, Center Freqency fc =' + string(fc_trial/(10^9)) + 'GHz, Zero Padding = '+ string(zero_padding) + 'N','FontName','Times')
%     end
        set(gca,'FontSize',10,'FontWeight','bold')
    colorbar
    
    velocities_20 = smoothdata(velocities_2);
    figure(2)
    subplot(1,2,d)
    plot(time_array,velocities_20,'b')
    % hold on
    % plot(time_array,velocities_3,'r')
    % hold off
    xlabel('\bf{Time (s)}','FontName','Times');
    xlim([-2 max(time_array)]+2)
    ylim([0 10])
    ylabel('\bf{Velocity (m/s)}','FontName','Times')
    
    %saveas(figure1,'OUTPUT\Signal.png')
    %saveas(figure2,'OUTPUT\TP_variation\EntireRawDopplerFreq_0.2.png')
    %saveas(figure5,'OUTPUT_CW_velocity\ZeroPadding_variation\Output_ZeroPadding' + string(zero_padding) + '.png')
end