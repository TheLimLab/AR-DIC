% Plot output of 'main_AR_DIC_process.m'
% 
% Run this script after main_AR_DIC_process.m to plot results.
% Run in sections (Ctrl + Enter) to view plots
% Basic plots:
%   Displacement magnitude map
%   Velocity vector field
%   Maximum and minimum displacement vs. time
%   Trace of velocity across frame through maximum velocity
%   Contraction volume vs. time
%   Contraction frequency heat map
%   Contraction magnitude heat map
%   Centroids of large contracting areas
%   Accumulative spatial contraction map
%   Smallest independent regions
%   Raw strain vs. time at maximum strain location
%   Principal strain vs. time at maximum strain location
%   Raw strain field
%   Principal strain field
%   Power spectrum of reference image
%   Synthetic topography image to export
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% MAIN_AR_DIC
% MAIN_AR_DIC_PROCESS

load('cmap_custom.mat')%load strain color map

%%
%------------------Displacement, velocity----------------------------------
contour_frame=plot_contour(stretch2,'disp',ncont,mapstyle); %displacement contour
quiver_frame=plot_quiver(stretch2,'vel',10,mapstyle,[0,20]); %velocity quiver

%supplementary:
% surf_frame=surfplot(stretch2,mapstyle);
% threshold_frame=plot_threshold_contour(stretch2,ncont,mapstyle);
% video=xy_quiver(stretch2,1);
% binary_video=plot_binary_mask(stretch2);
%--------------------------------------------------------------------------

%%
%----------------plot max displacement, bpm, velocity trace-----------------
figure %max displacement
hold on
plot_disp_location(stretch2,'max');
plot_disp_location(control2,'max');
hold off
title('Displacement at maximum');

figure %min displacement
hold on
plot_disp_location(stretch2,'min');
plot_disp_location(control2,'min');
hold off
title('Displacement at minimum');

figure %plot velocity trace across frame at maximum velocity location
hold on
svel_row=stretch2.frame_holder{1,stretch2.max_velocity_frame}.disp_max_index(1);
svel=stretch2.frame_holder{1,stretch2.max_velocity_frame}.vel_mat(svel_row,:);
cvel_row=control2.frame_holder{1,control2.max_velocity_frame}.disp_max_index(1);
cvel=control2.frame_holder{1,control2.max_velocity_frame}.vel_mat(svel_row,:);
plot(stretch2.x_mat(1,1:end),svel)
plot(control2.x_mat(1,1:end),cvel)
hold off
title(['Contraction velocity (\mum/s)']);
ylabel('velocity (\mum/s)')
xlabel('(\mum)');

%supplementary:
% figure
% hold on
% svel_row_min=stretch2.frame_holder{1,stretch2.min_velocity_frame}.disp_min_index(1);
% svel_min=stretch2.frame_holder{1,stretch2.min_velocity_frame}.vel_mat(svel_row_min,:);
% cvel_row=control2.frame_holder{1,control2.min_velocity_frame}.disp_min_index(1);
% cvel=control2.frame_holder{1,control2.max_velocity_frame}.vel_mat(svel_row_min,:);
% plot(svel_min)
% plot(cvel)
% hold off
%--------------------------------------------------------------------------

%%
%---------------------plot contraction volume------------------------------
cv_hand_contract=plot_contraction_volume(stretch2);
cv_hand_control=plot_contraction_volume(control2);

%%
%-----------------------Heat maps, Large area centroids--------------------
figure_handle=heatmap_frequency(stretch2,ncont,'jet');%frequency heat map
magnitude_heatmap=heatmap_magnitude(stretch2,ncont,'jet');%magnitude heat map
area_cent_hand=plot_centroid_areathreshold(stretch2,4840); %4840

%%
%-Accumulative spatial contraction (ASC) map, Smallest independent regions-
ASC_map(stretch2.tree,'3D',0);%ASC map
plot_node_array(leaves_s100,stretch2);%Smallest independent regions
%--------------------------------------------------------------------------

%%
%------------------------Strain trace plots--------------------------------
%raw strains at maximum location
plot_strain(strains_s2,'raw',xind_s2,yind_s2,stretch2.timestep);
plot_strain(strains_c2,'raw',xind_c2,yind_c2,control2.timestep);

%principal strains at maximum location
plot_strain(strains_principal_s2,'principal',xind_s2,yind_s2,stretch2.timestep);
plot_strain(strains_principal_c2,'principal',xind_c2,yind_c2,control2.timestep);

%supplementary
%plot_strain(strains_principal_c2,'principal',xind_min_c2,yind_min_c2,control2.timestep);
%plot maximum compression
%plot_strain(strains_principal_s2,'principal',xind_min_s2,yind_min_s2,stretch2.timestep);

%%
%--------------------------Raw Strain maps---------------------------------
%plot strain maps EXX EYY EXY
mapstyle2=cmap_comp_tens;% mapstyle2='hot';
frame_s=stretch2.max_disp_frame;
figure
contourf(stretch2.x_mat,stretch2.y_mat,strains_s2{1,frame_s}.Exx,ncont,'LineStyle','none');
title(['Strain ',char(949),'_x_x']);
axis ij
axis([0 stretch2.xdim 0 stretch2.ydim])
caxis([-0.1,0.1]) %caxis([-max_val_s2,max_val_s2])
colormap(mapstyle2)
colorbar;
figure
contourf(stretch2.x_mat,stretch2.y_mat,strains_s2{1,frame_s}.Eyy,ncont,'LineStyle','none');
title(['Strain ',char(949),'_y_y']);
axis ij
axis([0 stretch2.xdim 0 stretch2.ydim])
caxis([-0.1,0.1])
colormap(mapstyle2)
colorbar;
figure
contourf(stretch2.x_mat,stretch2.y_mat,strains_s2{1,frame_s}.Exy,ncont,'LineStyle','none');
title(['Strain ',char(949),'_x_y']);
axis ij
axis([0 stretch2.xdim 0 stretch2.ydim])
caxis([-0.1,0.1])
colormap(mapstyle2)
colorbar;
%--------------------------------------------------------------------------

%%
%--------------------Principal Strain maps---------------------------------
%strain map S1, S2
mapstyle2=cmap_comp_tens;% mapstyle2='hot';
frame_s=stretch2.max_disp_frame;
figure
contourf(stretch2.x_mat,stretch2.y_mat,strains_principal_s2{1,frame_s}.s1,ncont,'LineStyle','none');
title(['Strain ',char(949),'_1']);
axis ij
axis([0 stretch2.xdim 0 stretch2.ydim])
caxis([-0.1,0.1]) %caxis([-max_val_s2,max_val_s2])
colormap(mapstyle2)
colorbar;
figure
contourf(stretch2.x_mat,stretch2.y_mat,strains_principal_s2{1,frame_s}.s2,ncont,'LineStyle','none');
title(['Strain ',char(949),'_2']);
axis ij
axis([0 stretch2.xdim 0 stretch2.ydim])
caxis([-0.1,0.1])
colormap(mapstyle2)
colorbar;
figure
contourf(stretch2.x_mat,stretch2.y_mat,strains_principal_s2{1,frame_s}.max_eng_shear,ncont,'LineStyle','none');
title(['Strain ',char(949),'_shear']);
axis ij
axis([0 stretch2.xdim 0 stretch2.ydim])
caxis([-0.1,0.1])
colormap(mapstyle2)
colorbar;

%supplementary:
% figure
% for nn=1:stretch2.num_frame
% mapstyle2=cmap_comp_tens;% mapstyle2='hot';
% frame_s=nn;
% contourf(stretch2.x_mat,stretch2.y_mat,strains_s2{1,frame_s}.Exx,ncont,'LineStyle','none');
% title(['Strain ',char(949),'_x_x']);
% axis ij
% axis([0 stretch2.xdim 0 stretch2.ydim])
% caxis([-0.1,0.1]) %caxis([-max_val_s2,max_val_s2])
% colormap(mapstyle2)
% colorbar;
% strain_frame(nn)=getframe(gcf);
% end
%--------------------------------------------------------------------------

%%
%--------------------Plot reference power spectrum-------------------------
val_image=imread('Series80_ref_frame_stretch.tif'); %reference frame
% val_image=imread('Ref_frame_no_contraction_series015.tif');
caxis_lim=[0,8.5];
plot_power_spectrum(val_image,scale100,caxis_lim);
%--------------------------------------------------------------------------

%%
%------------------Intuitive topography map in Tangram---------------------
%input parameters
save_images=0;
bits=16; %image bitdepth to save
high_alpha=0.82;  %.82 (Mt. Diablo elevation)
low_alpha=.93;   %.93  (Sea level)
tile_size=256; %standard tile size of Tangram tile (256 x 256 px)
opt='all'; %'all' or 'maxmin'
row=5;
frame_holder={[51,100],[111,160]}; %{[51,100],[1,50]};
[RGBmat,Alpha]=export_to_map(data,frame_holder,row,high_alpha,low_alpha,tile_size,opt);
if save_images==1
    imwrite(RGBmat,'testrgb.png','png','Alpha',Alpha,'BitDepth',bits);
    imwrite(rgb2gray(RGBmat),'testgray.png','png','BitDepth',bits);
end
imshow(RGBmat)
%--------------------------------------------------------------------------

%%
%----------------------------FFT BPM---------------------------------------
% figure
% hold on
disp_max=plot_disp_location(stretch2,'max'); %O' 
% [~]=plot_disp_location(stretch2,[34 22]); %O''

% figure %plot disp at O', O'', O'''
% hold on
% % [~]=plot_disp_location(stretch2,'max'); %max disp
% [~]=plot_disp_location(stretch2,[22 33]); %O' RGB: 121,20,54
% [~]=plot_disp_location(stretch2,[34 22]); %O'' 2,38,60
% [~]=plot_disp_location(stretch2,[30 12]); %O''' 244,195,137
% hold off

% disp_max_control=plot_disp_location(control2,'max');
% [~]=plot_disp_location(stretch2,[34 24]);
% hold off

%--------FFT method for BPM--------------------
disp_max(disp_max<1)=0;
Fs=1/stretch2.timestep;% Sampling frequency
[Power1,freq]=fft_plot(disp_max,Fs);
[BPM,f_index]=bpm_from_fft(Power1,freq);
f1=freq(f_index);
%-----------------------------------------------

%-----------Peak counter method for BPM-----------
% [~,numpeak]=bwlabel(disp_max>0.14);% (disp_max<0.5)&(disp_max>0.14)
% BPM_peak_count=60*(numpeak/(length(disp_max)*timestep));
%-------------------------------------------------
%--------------------------------------------------------------------------

%%
% play_animation(contour_frame,1,3)
% play_animation(quiver_frame(83),1,3)
% play_animation(surf_frame,1,3)
% play_animation(threshold_frame,1,3)
% play_animation(binary_video,1,3)
% play_animation(surf_frame(48),1,1)%single frame

%%
%%%---------------To save video--------------------------------  
% mov=binary_video;
% movie2avi(mov, 'export_video')
%%%-------------------------------------------------------------

%EOF