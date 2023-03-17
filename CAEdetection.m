%---------------------------------------
% Automatic CAE detection algorithm
% Barbara Bica, Jorge Figueroa, Lisa Salgueiro, Rojan Aslani
% Superior Institute of Engineering of Porto (ISEP)
%---------------------------------------
clear all
close all
clc
clf

pkg load signal
pkg load io

tic = tic ();

spikes_excel = 'C:\Users\aslan\OneDrive - Instituto Superior de Engenharia do Porto\ISEP\21-2\PSFISIO\PRJ-GroupB1\Octave\20220617\spikes.xlsx';

%---------------------------------------
%---------------------------------------
% INPUTTING
%---------------------------------------
%---------------------------------------
fs = 256; # sampling frequency
disp("Uploading file...")
EEG = importdata('1564C_Sever.txt');
col7 = EEG(:,7); % Fz channel
N = length(col7);
t_ini = (0 : 1/fs : ((N-1)*1/fs));

% Fast Fourier Transformation (time -> frequency)
disp("FFT...")
NFFT = 2^nextpow2(N) ;
f = 0 : fs/(NFFT-1) : fs ; %(NFFT-1) is the number of intervals between 0 and NFFT
EEG_f = fft(col7,NFFT) * 2/N ;
t = (0 : 1/fs : ((NFFT-1)*1/fs));

%---------------------------------------
%---------------------------------------
% ANALYSIS
%---------------------------------------
%---------------------------------------

% ------ PRE-PROCESSING ------
disp("Applying Chebychev filters...")
n = 9;

% Chebychev II band pass filter - Noise reduction
fc_inf = 0;
fc_sup = 35;
EEG_cheby2_bandpass = cheby2_bandpass(NFFT, fs, n, fc_inf, fc_sup );
EEG_f_filteredCheby2 = EEG_f'.* EEG_cheby2_bandpass; % Apply filter to EEG signal

% ------ Chebychev II band pass filter - separate 4 brain waves ------
##disp("Delta")
delta_fc_inf = 0.5;
delta_fc_sup = 4;
delta_filter_cheby2_bandpass = cheby2_bandpass(NFFT, fs, n, delta_fc_inf, delta_fc_sup );
delta_EEG_f_filtradoCheby2 = EEG_f'.* delta_filter_cheby2_bandpass; % Apply filter to EEG signal

##disp("Theta")
theta_fc_inf = 3.5;
theta_fc_sup = 8.5;
theta_filter_cheby2_bandpass = cheby2_bandpass(NFFT, fs, n, theta_fc_inf, theta_fc_sup );
theta_EEG_f_filtradoCheby2 = EEG_f'.* theta_filter_cheby2_bandpass; % Apply filter to EEG signal

##disp("Alpha")
alpha_fc_inf = 7.5;
alpha_fc_sup = 12.5;
alpha_filter_cheby2_bandpass = cheby2_bandpass(NFFT, fs, n, alpha_fc_inf, alpha_fc_sup );
alpha_EEG_f_filtradoCheby2 = EEG_f'.* alpha_filter_cheby2_bandpass; % Apply filter to EEG signal

##disp("Beta")
beta_fc_inf = 11.5;
beta_fc_sup = 36;
beta_filter_cheby2_bandpass = cheby2_bandpass(NFFT, fs, n, beta_fc_inf, beta_fc_sup );
beta_EEG_f_filtradoCheby2 = EEG_f'.* beta_filter_cheby2_bandpass; % Apply filter to EEG signal

%---------------------------------------
%---------------------------------------
% PROCESSING - spike detection
%---------------------------------------
%---------------------------------------
% ------ IFFT -> changing to time domain ------
disp("IFFT...")

delta_EEG_t_filtradoCheby2 = ifft(delta_EEG_f_filtradoCheby2, NFFT)*(NFFT/2);

beta_EEG_t_filtradoCheby2 = ifft(beta_EEG_f_filtradoCheby2, NFFT)*(NFFT/2) ;

alpha_EEG_t_filtradoCheby2 = ifft(alpha_EEG_f_filtradoCheby2, NFFT)*(NFFT/2) ;

theta_EEG_t_filtradoCheby2 = ifft(theta_EEG_f_filtradoCheby2, NFFT)*(NFFT/2) ;

# --------- Taking zeros ---------------
disp("Zero removal...")
matrix = delta_EEG_t_filtradoCheby2;
minimum = 2;
n = 200;
delta_EEG_t_clean = check_zeros(matrix,minimum,n);
##delta_EEG_t_clean = delta_EEG_t_filtradoCheby2;

% ------ CROSS CORRELATION IN TIME DOMAIN ------
% Upload Spikes matrix 
disp("Uploading spikes matrix...")

[~,sheet_names] = xlsfinfo(spikes_excel); % ~ indica se o ficheiro é lido pela funcao xlsread e é descartado
for i=1:numel(sheet_names(:,1))
  #numel() retorna o número de elementos de sheet_name
  data = xlsread(spikes_excel,sheet_names{i});
  #plot(t_ini(1:length(seg)), seg)
  spikes_matrix(i,1:length(data)) = data;
end


disp("Cross correlating...")
rs = [];
lsm = length(spikes_matrix(:,1));

for i = 1:lsm
  segment = nonzeros(spikes_matrix(i,:));
  [r, lags] = xcorr( delta_EEG_t_clean, segment);
  ra = r( ceil(length(r)/2) : end );
  rs (i,:) = ra;
end

rs1 = rs(:,2:end);
r = mean(rs1);

% ------ THRESHOLDING ------
sizeofwindow = 5000;
rest = rem (length(delta_EEG_t_clean), sizeofwindow);

xs = [];
pos_peaks = [];
peaks_positions = [];
for i = (1 : (length(r)-rest+1)/sizeofwindow)
  rsmall = abs(r((i-1)*sizeofwindow+1 : i*sizeofwindow));
  threshold = mean(rsmall);
  maximum  = max(rsmall);
  x = (100 * threshold) / maximum;
  xs (i) = x;
  
  for j = (1:length(rsmall))
    blah = rsmall(j)*100/maximum;
    
    if (blah > x)
      pos_peaks(j+(i*sizeofwindow)) = blah;
      peaks_positions(end+1) = j+(i*sizeofwindow);
    end
  end
  
  xx = ones(1,length(rsmall))*x;
  figure (i)
  plot(lags(1:length(rsmall)), abs(rsmall)/maximum*100 ), axis tight
  hold on
  plot(lags(1:length(rsmall)), xx), axis tight
  hold off
  
end
  
if x < 25
    disp('Epileptic discharge detected!')
else
    disp('No epileptic discharge.')
end

% low pass - for threshold 
y_c = r/max(r);
for element = 1:length(r)
    if y_c(element) <= threshold/max(r)
        y_c(element) = 0;
    end
end
disp(y_c)
figure(2) 
plot(x,y_c)

%---------------------------------------
%---------------------------------------
              % VISUALIZATION
%---------------------------------------
%---------------------------------------

##plot( t (1:length(delta_EEG_t_clean)), delta_EEG_t_clean)
##index_list = [1,8,20,21,22,23];
##
##count = 1;
##
##for element = 1:length(index_list)
##    count = count + 1 ;
##    disp(count)
##    if mod(count,2) == 0
##        %disp(index_list(element))
##        %disp(index_list(element+1))
##        disp('loop')
##        hold on
##        i1 = index_list(element);
##        i2 = index_list(element+1);
##        plot([i1,i1],[-1,1])
##        plot([i2,i2],[-1,1])
##        idx=x>i1&x<i2;
##        H=area(x(idx),y1(idx));
##        set(H,'FaceColor',[1 0 0]);
##        hold off
##    end
##end

##set(0, "defaulttextfontsize", 10)  % title
##set(0, "defaultaxesfontsize", 7)  % axes labels

disp("EOF.")

toc = toc (tic);
disp("Time of execution: ")
disp(toc)

