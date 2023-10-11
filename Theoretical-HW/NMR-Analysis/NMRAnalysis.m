load nmrdata;
time = timeseconds;
Fs = (time(2,1) - time(1,1))^(-1);
T = 1/Fs;
% swap = Fs;
% Fs = T;
% T = swap;
L = 400001;

figure(1);
plot(time,EMF)
title('Figure 1: Raw signal');
xlabel('t (seconds)');
ylabel('EMF (mV)');

Y = fft(EMF);

Trans = abs(Y/L);

figure(2);
plot(Trans, '-r');
title('Figure 2:A plot of scaled FFT');
xlabel('Array Index');
ylabel('FFT/L');

L = L - 1; % to get rid of non-integer values
Spectra = Trans(1:(L/2+1));
Spectra(2:end-1)=Spectra(2:(end-1))*2;

f = (Fs*(0:(L/2))/L)';

figure(3);
plot(f,Spectra,'-r');
title('Figure 3: Single-Sided Amplitude Spectrum of NMR Data');
xlabel('Frequency (Hz)');
ylabel('|P(f)|');

figure(4);
plot(f(100000:107000),Spectra(100000:107000));
title(sprintf('Figure 4: Single-Sided Amplitude Spectrum of NMR Data \n zoomed in to view just the region around the peaks'));
xlabel('frequency (Hz)');
ylabel('|P(f)|');





