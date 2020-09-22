clear all
close all

  data = audioread('single_target_FMCW.wav');
% data = audioread('walking_opposite_ver3.wav');
% data = recordedData;

fs=44.1*10^3; % to read file
FMCW_range(data)

function [times, ranges] = FMCW_range(data)
    %Constants
    c = 3*10^8;   %m/s
    start_f = 2.408*10^9;   %Hz
    end_f = 2.495*10^9;   %Hz
    data = -data;
    sync_pulse = data(:,2);
    data_values = data(:,1);
    fs = 44100;

    %Realize the binary sync pulse %%ERROR: forget to (-1) the pulse
    for i=1:length(sync_pulse)
        if sync_pulse(i) > 0.1	%This threshold is used in order to avoid oscillations of the analog signal
           sync_pulse(i) = 1; 
        elseif sync_pulse(i) < - 0.1
            sync_pulse(i) = -1; 
        else
            sync_pulse(i)=0;
        end  
    end
%     upChirpStart = 0;   
%     upChirpEnd = 0;
%     for i=1:length(sync_pulse)-1
%         if sync_pulse(i + 1) == 1 && sync_pulse(i) <= 0
%             upChirpStart = i + 1 ;
%         elseif  sync_pulse(i) == 1 && sync_pulse(i + 1) <= 0
%             upChirpEnd = i;
%             if upChirpStart ~=  0 && upChirpEnd ~= 0
%                 diff = upChirpEnd - upChirpStart
%                 if diff < 700
%                     sync_pulse(upChirpStart:upChirpEnd) = 0 ;
%                 end
%                 upChirpStart = 0;
%                 upChirpEnd = 0;
%             end
%         end
%     end
    
    kd = 0 ; %% calculate number of down chirps

    % Calculate # of up-chirp
    k = 0;  %it indicates the # of up-chirp (# of rows of the matrix)
    for i = 2:length(data(:,1)) 
        if (sync_pulse(i-1) < 1 && sync_pulse(i) > 0)
            k = k+1;      
        elseif (sync_pulse(i-1) > -1 && sync_pulse(i) == -1)
            kd = kd+1;  
        end
    end

    % Matrix data
    N = 1100;
    updata_mat = zeros(k,N);% be careful
	upcolumn = 0;
    uprow = 1;
    
    downdata_mat = zeros(kd,N);
    downcolumn = 0;
    downrow = 1;
    
    % sort the data into two different matrices for up and down-chirp data
    for i = 1:(length(data(:,1))-1)
        if sync_pulse(i) == 1
            upcolumn = upcolumn + 1;
            updata_mat(uprow,upcolumn) = data_values(i);
            if sync_pulse(i+1) ~= 1
                uprow = uprow + 1;
                upcolumn = 0;
            end
        elseif sync_pulse(i) == -1
            downcolumn = downcolumn + 1;
            downdata_mat(downrow,downcolumn) = data_values(i);
            if sync_pulse(i+1) ~= -1
                downrow = downrow + 1;
                downcolumn = 0;
            end
        end
        
    end

    %Subtract the mean of each column from each column (MS clutter rejection)
    for i = 1:N
        b = sum(updata_mat(:,i))/length(updata_mat(:,i));
        updata_mat(:,i) = updata_mat(:,i) - b;
    end

    for i = 1:N
        b = sum(downdata_mat(:,i))/length(downdata_mat(:,i));
        downdata_mat(:,i) = downdata_mat(:,i) - b;
    end

    % Perform 2-pulse MTI
%     copy = data_mat;
%     for j = 2:k
%         for i = 1:N
%             data_mat(j,i) = copy(j,i)-copy(j-1,i); 
%         end
%     end

    % Perform 3-pulse MTI
    copy = updata_mat;
    for j = 3:k
        for i = 1:N
            updata_mat(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
        end
    end

    copy = downdata_mat;
    for j = 3:kd
        for i = 1:N
            downdata_mat(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
        end
    end
    %Perform IFFT
    for i = 1:length(updata_mat(:,1))
        upmat_time(i,:) = ifft(updata_mat(i,:),4*N);
    end
    for i = 1:length(updata_mat(:,1))
        downmat_time(i,:) = ifft(downdata_mat(i,:),4*N);
    end
    
    % take log of the range data 
    upmat_time = 20*log10(abs(upmat_time));
    downmat_time = 20*log10(abs(downmat_time));
    
    
    %Division of matrix by 2 and normalization
    upmat_time = upmat_time(:,1:(length(upmat_time(1,:))/2)); %divides the matrix by 2
    upmat_time = upmat_time - max(max(upmat_time));
    
    %Division of matrix by 2 and normalization
    downmat_time = downmat_time(:,1:(length(downmat_time(1,:))/2)); %divides the matrix by 2
    downmat_time = downmat_time - max(max(downmat_time));
    
    %calculate range resolution
    delta_f = end_f-start_f;
    delta_R = c/(2*delta_f); 

    %calculate max range
    % N_freq = (length(mat_time));
    R_max = N*delta_R/2;  
    [M,~] = size(upmat_time);
    
    T_max = length(data(:,1))/fs;
    range_array = linspace(0,R_max,N);
    time_array = linspace(0,T_max,M);

    % -- Velocity Calculation - Finite Different Method
    [~,I] = max(upmat_time');
    ranges = range_array(I)/2;
    %ranges = ranges(1:length(time_array));
    ranges_1 = smoothdata(ranges);
    times = time_array;
    for i=2:length(ranges_1)
        velocity(i)=(ranges_1(1,i)-ranges_1(1,i-1))/(time_array(i)-time_array(i-1));
    end
    velocity_1 = smooth(velocity);
    velocity_2 = smoothdata(velocity_1);
    
%     
% %    Plot signal
%     figure1 = figure;
%     plot(-data(:,1));
%     xlabel('Data sample number','FontName','Times');
%     xlim([0 4.5*10^4]);
%     ylim([-0.4 0.4]);
%     ylabel('Amplitude','FontName','Times');
%     title('Sampled down-converted data','FontName','Times');
%     set(gca,'FontSize',10,'FontWeight','bold');
% 
    figure2 = figure; hold on;
    plot(-data(:,2));
    plot(zeros(3*10^4,1),'--k','LineWidth',2);
    xlabel('Data sample number','FontName','Times');
    xlim([0 3*10^4]);
    ylim([-1.15 1.15]);
    ylabel('Amplitude','FontName','Times');
    title('Sync data','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    text1 = text(2.6*10^4,0.1,'Threshold','FontSize',10,'FontWeight','bold','FontName','Times');
    hold off;

    figure3 = figure; hold on;
    plot(sync_pulse);
    plot(zeros(3*10^4,1),'--k','LineWidth',2);
    plot(-data(:,1),'r')
    legend('Sync_Pulse',' ','Data')
    xlabel('Data sample number','FontName','Times');
    xlim([0 3*10^4]);
    ylim([-1.15 1.15]);
    ylabel('Amplitude','FontName','Times');
    title('Sync data','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    text2 = text(2.6*10^4,0.1,'Threshold','FontSize',10,'FontWeight','bold','FontName','Times');
    hold off;
 
    %%
    figure(5)
    subplot(1,2,1);
    imagesc(range_array,time_array,upmat_time,[-50,0])
    xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    xlim_inf = 0;
    xlim_sup = 100;
    xlim([xlim_inf xlim_sup])
    %text1 = text(70,37,'3-pulse MTI','FontSize',15,'FontWeight','bold','FontName','Times');
    colorbar
    % downchirp range data
    subplot(1,2,2);
    imagesc(range_array,time_array,downmat_time,[-50,0])
    xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    xlim_inf = 0;
    xlim_sup = 100;
    xlim([xlim_inf xlim_sup])
    %text1 = text(70,37,'3-pulse MTI','FontSize',15,'FontWeight','bold','FontName','Times');
    colorbar
    
    diff_mat = (downdata_mat(1:min(kd,k),:)-flip(updata_mat(1:min(kd,k),:)));
    diff_mat = 20*log10(abs(fft(diff_mat)));

        %Division of matrix by 2 and normalization
%     upmat_time = upmat_time(:,1:(length(upmat_time(1,:))/2)); %divides the matrix by 2
%     upmat_time = upmat_time - max(max(upmat_time));
    f=2.43*10^9;% 2.4 GHz
    v_max = c*fs/2/(2*f);

    velocity_array = linspace(0,v_max,length(diff_mat));
    max_value = max(max(diff_mat)); %max value of the matrix
    diff_mat =diff_mat-max_value; %normaliz
    figure(); imagesc(range_array,time_array,diff_mat,[-50,0])
    colorbar
    %% 
end
