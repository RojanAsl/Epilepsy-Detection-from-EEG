clear all
close all
clc
clf

pkg load io
pkg load signal

fs = 256; # sampling frequency
disp("Uploading file...")
EEG = importdata('1564C_Sever.txt');
col7 = EEG(:,7); % Fz channel
N = length(col7);

disp("Processing...")
NFFT = 2^nextpow2(N) ;
EEG_f = fft(col7 , NFFT) * 2/N ;

# Cheby II 18th order
n = 9;
delta_fc_inf = 0.5;
delta_fc_sup = 4;
delta_filter_cheby2_bandpass = cheby2_bandpass(NFFT, fs, n, delta_fc_inf, delta_fc_sup );

delta_EEG_f_filtradoCheby2 = EEG_f'.* delta_filter_cheby2_bandpass; % Apply filter to EEG signal
delta_EEG_t_filtradoCheby2 = ifft(delta_EEG_f_filtradoCheby2, NFFT)*(NFFT/2);

t = (0 : 1/fs : ((NFFT-1)*1/fs));

##figure (1)
##plot(t, delta_EEG_t_filtradoCheby2), axis tight
##title('Delta wave')

#segment = delta_EEG_t_filtradoCheby2 (floor(649.2*fs) : ceil(650.3*fs)); %1316C_Moderate
segment = delta_EEG_t_filtradoCheby2 (floor(372.7*fs) : ceil(373.5*fs)); % from 1564C_Sever

baseFileName = "C:\Users\aslan\OneDrive - Instituto Superior de Engenharia do Porto\ISEP\21-2\PSFISIO\PRJ-GroupB1\Octave\20220618\spikes.xlsx";
xlswrite(baseFileName, segment,'1564C_Sever_2');

#PLOT
##segment_ = xlsread(baseFileName, '1564C_Sever_2');
##figure (2)
##plot( t(1:length(segment_)), segment_), axis tight
##title('Segment')
