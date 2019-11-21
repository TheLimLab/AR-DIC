function plot_power_spectrum(input_image,scale,caxis_lim)
% PLOT_POWER_SPECTRUM Plots power spectrum of input_image.
% Power spectrum plot by taking the 2D FFT, shifting the zero component to
% the middle and plotting with a log10 transform.
%
% Usage: plot_power_spectrum(input_image,scale,caxis_lim)
% 
% Inputs: input_image is an nxmx3 array containing the image data.
% input_image may be obtained from function imread. scale (scalar) defines
% the pixel scale (for example pixels per micrometer). caxis_lim defines
% the color axis limits in the form [lower_limit upper_limit].
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln

figure
Fs_y=scale;
[ydim, xdim,~]=size(input_image); % number of pixels
dFx=scale/xdim;
dFy=Fs_y/ydim;
Fx=(-scale/2:dFx:scale/2-dFx)'; % cycles per micrometer
Fy=(-Fs_y/2:dFy:Fs_y/2-dFy)';
input_image=rgb2gray(input_image);
F_val=fft2(input_image);
shiftF=fftshift(F_val);
Log_scale=log10(1 + abs(shiftF));
image(Fx,Fy,Log_scale,'CDataMapping','scaled'); %x-y in units of 1/um (cycles per micrometer)
axis('square')
colormap('jet')
colorbar;
caxis(caxis_lim);
axis on;
xlabel('(\mum^-^1)')
ylabel('(\mum^-^1)')
title('log_1_0(|F|)');