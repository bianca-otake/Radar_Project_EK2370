%SAR script - Build your own radar
%Range migration algorithm (RMA)

%Initialize MATLAB
clear all
close all
clc

%Read file
[Y,Fs] = audioread('SAR_Test_File.m4a');

%%
%Constants definition
c = 3e8;    %m/s
fstart = 2.408*10^9;    %Hz
fstop = 2.495*10^9;
fc = 2.43*10^9;
BW = fstop - fstart;
Tsample = 20e-3;   %sec
Trp = 250e-3;  %minimum range profile time duration to acquire
                    %the range profile at one location
N = Fs*Tsample; %# of sample in each pulse
Nrp = Fs*Trp;   %# of sample between each range profile

%Seperate trigger and backscatter data
data = -1*Y(:,1);
trig = -1*Y(:,2);
%%

%Parse the data and sync signal for each range profile
rp_start = abs(trig) > mean(abs(trig));
count = 0;
% RP = zeros(24,11025);
% RP_trig= zeros(24,11025);
for i = Nrp+1 : (size(rp_start,1) - Nrp)
    if(rp_start(i) == 1 && sum(rp_start(i-Nrp : i-1)) == 0)
        count = count + 1;
        RP(count,:) = data(i : i+Nrp-1);
        RP_trig(count,:) = trig(i:i+Nrp-1);
    end
end
disp('done with parsing')

%% 
BW = (fstop - fstart)*1.2;
N = Fs*Tsample; %# of sample in each pulse
Nrp = Fs*Trp;   %# of sample between each range profile

thresh = 0.08;
clear sif; clear sif_h

for j = 1:size(RP,1)
    clear SIF;
    SIF = zeros(N,1);
    start = (RP_trig(j,:) > thresh);
    count = 0;
    for i = 12:(size(start,2)-2*N)   %%%%%
        [Y,I] = max(start(1,i:i+2*N));
        if (mean(start(i-10:i-2)) == 0 && I == 1)
            count = count + 1;
            SIF = RP(j,i:i+N-1)' + SIF;
        end
    end
    SI = SIF/count;
    FF = ifft(SI);
    clear SI;
    sif(j,:) = fft(FF(size(FF,1)/2 + 1:size(FF,1)));
end

%%

%Replace every NaN value with 1e-30
sif(isnan(sif)) = 1e-30;

%Other parameter definition
cr = BW/Tsample;                %chirp rate
Rs = 0;                         %down range distance to the scene center
lambda = c/fc;
delta_x = lambda/2;             %radar spacing bw each range profile aquisition data
L = delta_x*(size(sif,1));      %total SAR aperture length
Xa = linspace(-L/2,L/2,(L/delta_x));    %cross-ange radar position on L
time = linspace(0,Tsample,size(sif,2));
Kr = linspace((4*pi/c)*(fc-BW/2),(4*pi/c)*(fc+BW/2),size(time,2));

%MS clutter rejection
for j = 1:size(sif,1)
    sif(j,:) = sif(j,:) - mean(sif,1);
end

%Hanning window application to clean the spectrum
clear N;
N = size(sif,2);
H = zeros(1,N);
for i = 1:N
    H(i) = 0.5 + 0.5*cos(2*pi*(i-N/2)/N);
end

for i = 1:size(sif,1)
    sif_h(i,:) = sif(i,:).*H;
end
sif = sif_h;
%
%Phase of the SAR data matrix plot before along track FFT
fig_count=1;
figure(fig_count);
S_image= angle(sif);
imagesc(Kr, Xa, S_image);
colormap('default');
xlabel('K_r(rad/m)');
ylabel('SAR Position, Xa(m)');
colorbar;
fig_count = fig_count+1;

%%
zpad= 2048;
s_zeros= zeros(zpad, size(sif,2));
for i= 1:size(sif,2)
    index = round((zpad-size(sif,1))/2);
    s_zeros(index+1:(index + size(sif,1)),i) = sif(:,i);
end
sif = s_zeros;

%FFT along 1 dimension
S = fftshift(fft(sif,[], 1), 1);

%Define Kx variation
Kx= linspace((-pi/delta_x), (pi/delta_x), (size(S,1)));

%Magnitude after along track FFT
figure(fig_count);
S_image= 20*log10(abs(S));
imagesc(Kr, Kx, S_image, [max(max(S_image))-40,max(max(S_image))]);
colormap('default');
xlabel('K_r(rad/m)');
ylabel('K_x(rad/m)');
colorbar;
fig_count= fig_count+ 1;

%Phase after along track FFT'
clear S_image
figure(fig_count);
S_image= angle(S);
imagesc(Kr, Kx, S_image);
colormap('default');
xlabel('K_r(rad/m)');
ylabel('K_x(rad/m)');
colorbar;
fig_count= fig_count+ 1;

%%
Rs = 0;
S_matched = S;

kstart=floor(min(Kr)-20);
kstop=ceil(max(Kr));
Ky_e=linspace(kstart,kstop,1024);
count = 0;
for ii = 1:zpad;
    count = count + 1;
    Ky(count,:) = sqrt(Kr.^2 -Kx(ii)^2);
    S_st(count,:) = (interp1(Ky(count,:), S_matched(ii,:), Ky_e));
end

S_st(find(isnan(S_st))) = 1E-30;


figure(fig_count);
S_image= angle(S_st);
imagesc(Ky_e, Kx, S_image);
colormap('default');
xlabel('K_y(rad/m)');
ylabel('K_x(rad/m)');
colorbar;
fig_count= fig_count+ 1;
%%
v = ifft2(S_st,(size(S_st,1)*4),(size(S_st,2)*4));
bw= c*(kstop-kstart)/(4*pi);
max_range= (c*size(S_st,2)/(2*bw));
figure(fig_count);
S_image= v;
S_image= fliplr(rot90(S_image));
cr1 = -25; % depends on the Kx of the StoltInterpolation
cr2 = 25; % depends on the Kx of the StoltInterpolation
dr1 = 1;
dr2 = 100;

dr_index1 = round((dr1/max_range)*size(S_image,1));
dr_index2 = round((dr2/max_range)*size(S_image,1));
cr_index1 = round(( (cr1+zpad*delta_x/(2*1)) / (zpad*delta_x/1))*size(S_image,2));


cr_index2 = round(((cr2+zpad*delta_x/(2*1)) / (zpad*delta_x/1))*size(S_image,2));
trunc_image= S_image(dr_index1:dr_index2,cr_index1:cr_index2);
downrange = linspace(-1*dr1,-1*dr2, size(trunc_image,1));
crossrange= linspace(cr1, cr2, size(trunc_image, 2));
for ii = 1:size(trunc_image,2)
trunc_image(:,ii) = (trunc_image(:,ii)').*(abs(downrange*1)).^(3/2);
end
trunc_image= 20 * log10(abs(trunc_image));
imagesc(crossrange, downrange, trunc_image, [max(max(trunc_image))-40, max(max(trunc_image))-0]);
colormap('default');
axis equal;
colorbar;


