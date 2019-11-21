function [BPM,f_index]=bpm_from_fft(P1,f)
% BPM_FROM_FFT Calculates BPM from output of fft_plot. 
%
% Usage: [BPM,f_index]=bpm_from_fft(P1,f)
%
% Inputs: P1, a 1xn array from the power spectrum of fft_plot. f, a 1xn
% array containing the frequencies from fft_plot
% 
% Returns: vector of sorted BPMs in descending order, BPM, and indices to
% original frequency, f_index.
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% FFT_PLOT

[~,f_index]=sort(P1,'descend');
Hz=f(f_index);
BPM=Hz*60;%convert Hz to BPM