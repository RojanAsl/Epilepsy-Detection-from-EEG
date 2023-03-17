function filtroCheby2 = cheby2_bandpass(NFFT, fs, n, fc_inf, fc_sup)
% NFFT
% fs = sample frequency
% n = filter order
% Rs = Ripple
% fc_inf = low cut frequency
% fc_sup = high cut frequency

w=0:2*pi/(NFFT-1):2*pi;
z=exp(1i*w);

f=0:fs/(NFFT-1):fs;     %(NFFT-1) é o número de intervalos entre 0 e NFFT

Rs = 50;
num=1;
den=1;

[z_ch,p_ch,k_ch] = cheby2( n , Rs ,([fc_inf fc_sup]./(fs./2)));


for i=1:length(z_ch)
  num = (z - z_ch(i)) .* num; %numerador
  den = (z - p_ch(i)) .* den; %denominador
end

filtroCheby2 = abs (k_ch .* num ./ den); % Create the filter

end
