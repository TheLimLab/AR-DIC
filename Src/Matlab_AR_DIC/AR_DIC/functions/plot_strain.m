function plot_strain(frames,type,row_ind,col_ind,step)
% PLOT_STRAIN Plots strain vs. time at specified row and column index.
%
% Usage: plot_strain(frames,type,row_ind,col_ind,step)
%
% When type is 'principal', the function plots the maximum strain, minimum
% strain, and maximum shear strain. When type is 'raw' the function plots
% strain in the horizontal (Exx), vertical (Eyy), and shear (Exy)
% directions.
% 
% Inputs: frames a 1xn cell array of strain structures, type (string)
% either 'principal' or 'raw' to plot principal or raw strains
% respectively. row_ind (scalar) row index of frame to plot, col_ind
% (scalar) column index of frame to plot.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also: 
% VIDOBJ.CALC_STRAIN
% VIDOBJ/STRAIN_TO_CELL
% VIDOBJ/TRANSFORM_STRAIN

%initialize vectors for holding strain values
num_frame=length(frames);
primary_strain=zeros(1,num_frame);
secondary_strain=primary_strain;
shear_strain=primary_strain;
t=(1:num_frame).*step;%time vector

switch type
    case 'principal'
        primary_angle=primary_strain;
        for n=1:num_frame
            %get strain values from all frames
            primary_strain(n)=frames{n}.s1(row_ind,col_ind);
            secondary_strain(n)=frames{n}.s2(row_ind,col_ind);
            shear_strain(n)=frames{n}.max_eng_shear(row_ind,col_ind);
            primary_angle(n)=frames{n}.theta(row_ind,col_ind);
        end
        figure %plot principal angle
        plot(t,primary_angle)
        title('Principal angle')
        xlabel('Time (s)')
        ylabel('Principal angle (degrees)')
        principal_title='Principal strains';
        primary_label=[char(949),'_1'];
        secondary_label=[char(949),'_2'];
        shear_label=[char(949),'_s_h_e_a_r _m_a_x'];
    case 'raw'
        for n=1:num_frame
            %collect raw strain values in X, Y, and shear directions
            primary_strain(n)=frames{n}.Exx(row_ind,col_ind);
            secondary_strain(n)=frames{n}.Eyy(row_ind,col_ind);
            shear_strain(n)=frames{n}.Exy(row_ind,col_ind);
        end
        principal_title='Raw strains';
        primary_label=[char(949),'_x_x'];
        secondary_label=[char(949),'_y_y'];
        shear_label=[char(949),'_x_y'];
end

%plot 3 subfigures for primary strain directions and shear
figure
hold on
subplot(3,1,1)
plot(t,primary_strain)%E1
ylabel(primary_label)
title(principal_title)

subplot(3,1,2)
plot(t,secondary_strain)%E2
ylabel(secondary_label)

subplot(3,1,3)
plot(t,shear_strain)%Eshear
ylabel(shear_label)
xlabel('Time (s)')
hold off