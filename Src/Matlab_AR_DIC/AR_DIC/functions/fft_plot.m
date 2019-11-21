function [P1,f]=fft_plot(in_vect,Fs)
% FFT_PLOT Takes FFT of in_vect and plots FFT vs frequency.
%
% Usage: [P1,f]=fft_plot(in_vect,Fs)
%
% Input: 1xn array in_vect, Fs (scalar) sampling frequency.
%
% Returns: vectors power spectrum, P1, and frequency, f.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% BPM_FROM_FFT

L=length(in_vect);% Length of signal
Y=fft(in_vect);
P2=abs(Y/L);
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1); %remap power spectrum
f=Fs*(0:(L/2))/L;
figure
plot(f,P1)
xlabel('Frequency (Hz)')
ylabel('Power')