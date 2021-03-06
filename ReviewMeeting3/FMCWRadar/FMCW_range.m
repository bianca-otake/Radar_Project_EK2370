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
        if(sync_pulse(i) > 0.1)	%This threshold is used in order to avoid oscillations of the analog signal
           sync_pulse(i) = 1; 
        else
            sync_pulse(i)=0;
        end  
    end

    % Calculate # of up-chirp
    k = 0;  %it indicates the # of up-chirp (# of rows of the matrix)
    for i = 2:length(data(:,1)) 
        if (sync_pulse(i-1) < 1 && sync_pulse(i) > 0)
            k = k+1;      
        end
    end

    % Matrix data
    N = 1100;
    data_mat = zeros(k,N);% be careful

    column = 0;
    row = 1;
    for i = 1:(length(data(:,1))-1) 
        if (sync_pulse(i) > 0)
            column = column + 1;
            data_mat(row,column) = data_values(i);
            if (sync_pulse(i+1) < 1)
                row = row + 1;
                column = 0;
            end
        end 
    end
    upChirpStart = 0;
    upChirpEnd = 0;
    for i=1:length(sync_pulse)-1
        if sync_pulse(i + 1) == 1 && sync_pulse(i) <= 0
            upChirpStart = i + 1 ;
        elseif  sync_pulse(i) == 1 && sync_pulse(i + 1) <= 0
            upChirpEnd = i;
            if upChirpStart ~=  0 && upChirpEnd ~= 0
                diff = upChirpEnd - upChirpStart
                if diff < 700
                    sync_pulse(upChirpStart:upChirpEnd) = 0 ;
                end
                upChirpStart = 0;
                upChirpEnd = 0;
            end
        end
    end
    %Subtract the mean of each column from each column (MS clutter rejection)
    for i = 1:N
        b = sum(data_mat(:,i))/length(data_mat(:,i));
        data_mat(:,i) = data_mat(:,i) - b;
    end

    % Perform 2-pulse MTI
%     copy = data_mat;
%     for j = 2:k
%         for i = 1:N
%             data_mat(j,i) = copy(j,i)-copy(j-1,i); 
%         end
%     end

    % Perform 3-pulse MTI
    copy = data_mat;
    for j = 3:k
        for i = 1:N
            data_mat(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
        end
    end

    %Perform IFFT
    for i = 1:length(data_mat(:,1))
        mat_time(i,:) = ifft(data_mat(i,:),4*N);
    end
    mat_time = 20*log10(abs(mat_time));

    %Division of matrix by 2 and normalization
    mat_time = mat_time(:,1:(length(mat_time(1,:))/2)); %divides the matrix by 2
    mat_time = mat_time - max(max(mat_time));

    %calculate range resolution
    delta_f = end_f-start_f;
    delta_R = c/(2*delta_f); 

    %calculate max range
    % N_freq = (length(mat_time));
    R_max = N*delta_R/2;  

    T_max = length(data(:,1))/fs;
    range_array = linspace(0,R_max,2*N);
    time_array = linspace(0,T_max,k);

    %Plot signal
%     figure1 = figure;
%     plot(-data(:,1));
%     xlabel('Data sample number','FontName','Times');
%     xlim([0 4.5*10^4]);
%     ylim([-0.4 0.4]);
%     ylabel('Amplitude','FontName','Times');
%     title('Sampled down-converted data','FontName','Times');
%     set(gca,'FontSize',10,'FontWeight','bold');

%     figure2 = figure; hold on;
%     plot(-data(:,2));
%     plot(zeros(3*10^4,1),'--k','LineWidth',2);
%     xlabel('Data sample number','FontName','Times');
%     xlim([0 3*10^4]);
%     ylim([-1.15 1.15]);
%     ylabel('Amplitude','FontName','Times');
%     title('Sync data','FontName','Times');
%     set(gca,'FontSize',10,'FontWeight','bold');
%     text1 = text(2.6*10^4,0.1,'Threshold','FontSize',10,'FontWeight','bold','FontName','Times');
%     hold off;

%     figure3 = figure; hold on;
%     plot(sync_pulse);
%     plot(zeros(3*10^4,1),'--k','LineWidth',2);
%     xlabel('Data sample number','FontName','Times');
%     xlim([0 3*10^4]);
%     ylim([-1.15 1.15]);
%     ylabel('Amplitude','FontName','Times');
%     title('Sync data','FontName','Times');
%     set(gca,'FontSize',10,'FontWeight','bold');
%     text2 = text(2.6*10^4,0.1,'Threshold','FontSize',10,'FontWeight','bold','FontName','Times');
%     hold off;
    
    [~,I] = max(mat_time');
    ranges = range_array(I);
    ranges = ranges(1:length(time_array));
    times = time_array;

%     figure4 = figure;
%     imagesc(range_array,time_array,mat_time,[-50,0])
%     xlabel('Range(m)','FontName','Times')
%     ylabel('Time(s)','FontName','Times')
%     title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
%     set(gca,'FontSize',10,'FontWeight','bold');
%     xlim_inf = 0;
%     xlim_sup = 100;
%     xlim([xlim_inf xlim_sup])
%     %text1 = text(70,37,'3-pulse MTI','FontSize',15,'FontWeight','bold','FontName','Times');
%     colorbar
return


% saveas(figure1,'OUTPUT_FMCW_range\signal.png');
% saveas(figure2,'OUTPUT_FMCW_range\sync_signal.png');
% saveas(figure3,'OUTPUT_FMCW_range\sync_binary_signal.png');
% saveas(figure4,'OUTPUT_FMCW_range\ZeroPadding\output6.png');