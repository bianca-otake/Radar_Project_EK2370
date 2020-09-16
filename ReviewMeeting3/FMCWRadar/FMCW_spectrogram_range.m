function c = FMCW_spectrogram_range(data)
    %Constants
    c = 3*10^8;   %m/s
    start_f = 2.408*10^9;   %Hz
    end_f = 2.495*10^9;   %Hz
    data = -data;
    fs = 44100;
    sync_pulse = data(:,2);
    data_values = data(:,1);


    %Realize the binary sync pulse %%ERROR: forget to (-1) the pulse
    for i=1:length(sync_pulse)
        if(sync_pulse(i) > 0.1) 	%This threshold is used in order to avoid oscillations of the analog signal
           sync_pulse(i) = 1; 
        elseif sync_pulse(i) < -0.1
            sync_pulse(i) = -1;
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
    
%     % Calculate # of down-chirp
%     kd = 0;  %it indicates the # of down-chirp (# of rows of the matrix)
%     for i = 2:length(data(:,1)) 
%         if (sync_pulse(i-1) > -1 && sync_pulse(i) < 0)
%             kd = kd+1;      
%         end
%     end
%     
    % Matrix data
    N = 916;
    data_mat = zeros(k,N);% be careful
%     down_data_mat = zeros(kd,N);% be careful

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
%     column = 0;
%     row = 1;
%     for i = 1:(length(data(:,1))-1) 
%         if (sync_pulse(i) < 0)
%             column = column + 1;
%             down_data_mat(row,column) = data_values(i);
%             if (sync_pulse(i+1) > -1)
%                 row = row + 1;
%                 column = 0;
%             end
%         end 
%     end
    %Subtract the mean of each column from each column (MS clutter rejection)
    for i = 1:N
        b = sum(data_mat(:,i))/length(data_mat(:,i));
        data_mat(:,i) = data_mat(:,i) - b;
    end
%     for i = 1:N
%         b = sum(down_data_mat(:,i))/length(down_data_mat(:,i));
%         down_data_mat(:,i) = down_data_mat(:,i) - b;
%     end
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
%     % Perform 3-pulse MTI
%     copy = down_data_mat;
%     for j = 3:kd
%         for i = 1:N
%             down_data_mat(j,i) = copy(j,i) - copy(j-1,i) - (copy(j-1,i) - copy(j-2,i));
%         end
%     end
    %Perform IFFT
    for i = 1:length(data_mat(:,1))
        mat_time(i,:) = ifft(data_mat(i,:),4*N);
    end
    mat_time = 20*log10(abs(mat_time));

    %Division of matrix by 2 and normalization
    mat_time = mat_time(:,1:(length(mat_time(1,:))/2)); %divides the matrix by 2
    mat_time = mat_time - max(max(mat_time));

    %Perform IFFT
%     for i = 1:length(down_data_mat(:,1))
%         down_mat_time(i,:) = ifft(down_data_mat(i,:),4*N);
%     end
%     down_mat_time = 20*log10(abs(down_mat_time));
% 
%     %Division of matrix by 2 and normalization
%     down_mat_time = down_mat_time(:,1:(length(down_mat_time(1,:))/2)); %divides the matrix by 2
%     down_mat_time = down_mat_time - max(max(down_mat_time));

    
    %calculate range resolution
    delta_f = end_f-start_f;
    delta_R = c/(2*delta_f); 

    %calculate max range
    % N_freq = (length(mat_time));
    R_max = N*delta_R/2;  

    T_max = length(data(:,1))/fs;
    range_array = linspace(0,R_max,N*2);
    time_array = linspace(0,T_max,k);

    %Plot signal
    figure1 = figure;
    plot(-data(:,1));
    xlabel('Data sample number','FontName','Times');
    xlim([0 4.5*10^4]);
    ylim([-0.4 0.4]);
    ylabel('Amplitude','FontName','Times');
    title('Sampled down-converted data','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');

    hold on;
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
    xlabel('Data sample number','FontName','Times');
    xlim([0 3*10^4]);
    ylim([-1.15 1.15]);
    ylabel('Amplitude','FontName','Times');
    title('Sync data','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    text2 = text(2.6*10^4,0.1,'Threshold','FontSize',10,'FontWeight','bold','FontName','Times');
    hold off;
    Nmax = 4
    [Tp,~] = islocalmax(mat_time,2,'MaxNumExtrema',Nmax);
    imshow(Tp)
    ranges = zeros(size(Tp,1),Nmax);
    for i = 1:size(Tp,1)
       ranges(i,:) = range_array(Tp(i,:)) ;
    end
    figure()
    hold on
    for i = 1:Nmax 
        scatter(time_array,ranges(:,i),'o')
    end
     xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    ylim_inf = 0;
    ylim_sup = 100;
    ylim([ylim_inf ylim_sup])
    
    
    [~,I] = max(mat_time')
    range = range_array(I);
    range = range(range >1 );
    times = time_array(range >1 );
    
    figure();
    subplot(2,1,1);
    plot(time_array,range)

    figure();
    imagesc(range_array,time_array,mat_time,[-50,0])
    xlabel('Range(m)','FontName','Times')
    ylabel('Time(s)','FontName','Times')
    title('RTI with clutter rejection, f_{start} = 2.408 GHz , f_{stop} = 2.495 GHz','FontName','Times');
    set(gca,'FontSize',10,'FontWeight','bold');
    xlim_inf = 0;
    xlim_sup = 100;
    xlim([xlim_inf xlim_sup])
    colorbar

end

