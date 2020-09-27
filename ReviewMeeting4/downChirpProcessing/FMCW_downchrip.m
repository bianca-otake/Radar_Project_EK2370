clear all
close all

data = audioread('single_target_FMCW.wav');
% load('fmcw_velocity_accel.mat');
% data = r_array;
%   data = audioread('Range_Test_File.m4a');
% data = audioread('walking_opposite_ver3.wav');
% data = recordedData;

fs=44.1*10^3; % to read file
FMCW_range2(data)

function [times, ranges] = FMCW_range2(data)
    %Constants
    c = 3*10^8;   %m/s
    start_f = 2.408*10^9;   %Hz
    end_f = 2.495*10^9;   %Hz
    center_f = (start_f+end_f)/2; % center frequency
    lamdha = c/center_f; % waavelength
    data = -data; % inverse operation after sound card
    sync_pulse = data(:,2); % Syn pulse 
    data_values = data(:,1); % main data
    fs = 44100; % sampling frequency

    %Realize the binary sync pulse %%ERROR: forget to (-1) the pulse
    for i=1:length(sync_pulse)
        if(sync_pulse(i) > 0.1)	%This threshold is used in order to avoid oscillations of the analog signal
           sync_pulse(i) = 1; 
        
        else
            sync_pulse(i) = -1;
        end  
    end
    
    % Calculate # of up-chirp
    k = 0;  %it indicates the of up-chirp (# of rows of the matrix)
    m = 0;  %it indicates the of down-chirp (# of rows of the matrix)
    for i = 2:length(data(:,1)) 
        if (sync_pulse(i-1) < 1 && sync_pulse(i) > 0)
            k = k+1;
        elseif (sync_pulse(i-1) >= 1 && sync_pulse(i) < 0)
            m=m+1;
        else
        end
    end
    % Matrix data
   N = 882;%916; % Number of samples in during sync pulse ON
   data_mat = zeros(k+m,N);% be careful
  
   column = 0;
    row = 1;
    for i = 1:(length(data(:,1))-1)
        column = column + 1;
        if (sync_pulse(i) == 1) % for Up chirp data
            data_mat(row,column) = data_values(i); % data into matrix form
            
            if (sync_pulse(i+1) ~= 1)
                row = row + 1; % change to next row
                column = 0; % reinitializing coloumb to be zero
            end
            
        elseif (sync_pulse(i) == -1)% similary for down chirp data
            data_mat(row,column) = data_values(i); 
            if (sync_pulse(i+1) ~= -1 )
                row = row + 1;
                column = 0;
            end
        end
    end
    
    for i=1:length(data_mat)/2
        data_down(i,:) = data_mat(2*i-1,:); % Down chirp data
        data_up(i,:) = data_mat(2*i,:); % Up chirp data
        data_doppler(i,:) = (data_down(i,:)-data_up(i,:));
    end
    
    %Subtract the mean of each column from each column (MS clutter rejection)
    for i = 1:N
        b = sum(data_doppler(:,i))/length(data_doppler(:,i));
        data_doppler(:,i) = data_doppler(:,i) - b;
    end
    for i = 1:N
        b = sum(data_up(:,i))/length(data_up(:,i));
        data_up(:,i) = data_up(:,i) - b;
    end
    for i = 1:N
        b = sum(data_down(:,i))/length(data_down(:,i));
        data_down(:,i) = data_down(:,i) - b;
    end

    % Perform 2-pulse MTI
%     copy = data_mat;
%     for j = 2:k
%         for i = 1:N
%             data_mat(j,i) = copy(j,i)-copy(j-1,i); 
%         end
%     end

    % Perform 3-pulse MTI
    copy = data_doppler;
    for j = 3:k
        for i = 1:N
            data_doppler(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
        end
    end
    
    copy = data_up;
    for j = 3:k
        for i = 1:N
            data_up(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
        end
    end
    
    copy = data_down;
    for j = 3:k
        for i = 1:N
            data_down(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
        end
    end
    
    %Perform IFFT
    for i = 1:length(data_doppler(:,1))
        mat_time_dopper(i,:) = ifft(data_doppler(i,:),4*N);
    end
    for i = 1:length(data_up(:,1))
        mat_time_up(i,:) = ifft(data_up(i,:),4*N);
    end
    for i = 1:length(data_down(:,1))
        mat_time_down(i,:) = ifft(data_down(i,:),4*N);
    end
    
    mat_time_dopper = 20*log10(abs(mat_time_dopper));
    mat_time_up = 20*log10(abs(mat_time_up));
    mat_time_down = 20*log10(abs(mat_time_down));

    %Division of matrix by 2 and normalization
    mat_time_dopper = mat_time_dopper(:,1:(length(mat_time_dopper(1,:))/2)); %divides the matrix by 2
    mat_time_dopper = mat_time_dopper - max(max(mat_time_dopper));
    
    mat_time_up = mat_time_up(:,1:(length(mat_time_up(1,:))/2)); %divides the matrix by 2
    mat_time_up = mat_time_up - max(max(mat_time_up));
    
    mat_time_down = mat_time_down(:,1:(length(mat_time_down(1,:))/2)); %divides the matrix by 2
    mat_time_down = mat_time_down - max(max(mat_time_down));
    
    %[m,n]=size(mat_time_dopper);
    [m3,n3]= size(mat_time_down);
    for i=2:m3
        [~,n] = max(mat_time_down(i,:));
        [~,n1] = max(mat_time_up(i,:));
        if(abs(n-n1)>5) % max doppler frequency change can not be more that 5*fs/n3
            %n-n1=10;
            index(i)=index(i-1);
        else 
            index(i)=n-n1;
        end
        doppler_frequency(i)= (n-n1)*fs/n3; % mapping the index value to frequency while considering the max doppler frequency tobe fs/2
    end
    
    %calculate range resolution
    delta_f = end_f-start_f;
    delta_R = c/(2*delta_f); 

    %calculate max range
    R_max = N*delta_R/2; 
    velocity_array = lamdha*doppler_frequency/4; % same as used in CW RADAR
    n = length( velocity_array );
    for i=2:n
        if (abs(velocity_array(i)-velocity_array(i-1))>6) % Condition max velocity can not be more than 6 m/s for normal human
            velocity_array(i)=velocity_array(i-1);
        end
    end
    
    [M,~] = size(mat_time_dopper);
    T_max = length(data(:,1))/fs;
    range_array = linspace(0,R_max,N);
    time_array = linspace(0,T_max,M);
    
    % -- Velocity Calculation - Finite Different Method
    [~,I] = max(mat_time_up');
    [~,I2] = max(mat_time_dopper');
    ranges_up= range_array(I)/2;
    ranges_down = range_array(I2)/2;
    ranges_up_2 = smoothdata(ranges_up);
    ranges_down_2 = smoothdata(ranges_down);

    velocity = velocity_array;
    velocity_2 = smoothdata(velocity);
    velocity_3 = smooth(velocity_2);
    velocity_4 = smoothdata(velocity_3);
    velocity_5 = smooth(velocity_4);
    velocity_6 = smoothdata(velocity_5);

%%% Finite difference method    
%    for i=2:length(ranges_up_2)
%         velocity(i)=(ranges_up_2(1,i)-ranges_up_2(1,i-1))/(time_array(i)-time_array(i-1));
%     end
%     velocity_1 = smooth(velocity);
%     velocity_2 = smoothdata(velocity_1);
    
    
    figure(5)
    imagesc(range_array,time_array,mat_time_up,[-50,0]) % Up chirp plot
    %imagesc(range_array,time_array,mat_time_down,[-50,0])
    xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    xlim_inf = 0;
    xlim_sup = 100;
    xlim([xlim_inf xlim_sup])
    colorbar
     
    figure(6)
    %imagesc(range_array,time_array,mat_time_up,[-50,0])
    imagesc(range_array,time_array,mat_time_down,[-50,0]) % Down chirp plot
    xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    xlim_inf = 0;
    xlim_sup = 100;
    xlim([xlim_inf xlim_sup])
    colorbar
    
    figure(7)
    load('ref_data_FDM.mat')
    plot(time_array,zeros(size(velocity_5)));
    hold on
    plot(time_array,velocity_3);
    plot(time_array,Ref_velocity_FDM,'--','LineWidth',4);
    ylabel('Velocity (m/s)','FontName','Times');
    xlim([0 max(time_array)]);
    xlabel('Time(s)','FontName','Times')
    set(gca,'FontSize',10,'FontWeight','bold');
    legend('Zero crossing','Up-down chirp data','Finite Difference Method')
    grid on
    hold off
    
%     figure(10)
%     plot(time_array,zeros(size(velocity_5)));
%     hold on
%     plot(time_array,velocity);
%     ylabel('Velocity(m/s)','FontName','Times');
%     xlim([0 max(time_array)]);
%     xlabel('Time(s)','FontName','Times')
%     hold off
%  
    figure(8)
    plot(time_array,ranges_up_2);
    hold on
    plot(time_array,ranges_down_2);
    ylabel('Range(m)','FontName','Times');
    xlim([0 max(time_array)]);
    xlabel('Time(s)','FontName','Times')
    legend('Up-chirp data','Down-chirp data')
    set(gca,'FontSize',10,'FontWeight','bold');
    hold off
    grid on
%     
%     figure(9)
%     yyaxis left
%     plot(time_array,ranges_1);
%     ylabel('Range(m)','FontName','Times');
%     xlim([0 max(time_array)]);
%     xlabel('Time(s)','FontName','Times')
%     yyaxis right
%     plot(time_array,velocity_2);
%     ylabel('Velocity(m/s)','FontName','Times');
%     xlim([0 max(time_array)]);
%     xlabel('Time(s)','FontName','Times')
    



end


% saveas(figure1,'OUTPUT_FMCW_range\signal.png');
% saveas(figure2,'OUTPUT_FMCW_range\sync_signal.png');
% saveas(figure3,'OUTPUT_FMCW_range\sync_binary_signal.png');
% saveas(figure4,'OUTPUT_FMCW_range\ZeroPadding\output6.png');
