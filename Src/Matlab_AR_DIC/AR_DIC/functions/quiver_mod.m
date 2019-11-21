function q=quiver_mod(x_mat,y_mat,u_mat,v_mat,qscale,mapstyle,mincaxis,maxcaxis)
% QUIVER_MOD Modifies quiver plot with color coding of arrows. 
% Output is same as quiver but with scaled color data.
%
% Usage: q=quiver_mod(x_mat,y_mat,u_mat,v_mat,qscale,mapstyle,mincaxis,maxcaxis)
%
% Input: x_mat,y_mat,u_mat,v_mat,qscale are same as in quiver function.
% Briefly, these are the matrices of x values, y values, x magnitudes, y
% magnitudes, and qscale (scalar) which is the quiver arrow scale. mapstyle
% specifies the color map to use. Any of the Matlab color maps may be
% specified. mincaxis (scalar) defines the lower value of the colormap.
% maxcaxis (scalar) defines the maximum value of the color axis.
%
% Returns: quivergroup handle.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% QUIVER

q=quiver(x_mat,y_mat,u_mat,v_mat,qscale); %plot quiver and then modify with color
axis ij %match orientation of other plots
colormap(mapstyle)
%color change code adapted from http://stackoverflow.com/questions/29632430/quiver3-arrow-color-corresponding-to-magnitude
mags=sqrt(sum(cat(2,q.UData(:),q.VData(:)).^2,2)); %compute vector mag
mags=[mincaxis; mags ;maxcaxis]; %concatenate with min and max values to get proper color scale
currentColormap=colormap(gca);% get colormap
[~,~,ind]=histcounts(mags,size(currentColormap,1)); %get color for arrow
cmap=uint8(ind2rgb(ind(:),currentColormap)*255); %get RGB color
cmap(1,:,:)=[]; %remove extra data from scaling min
cmap(end,:,:)=[];%remove extra data from scaling max
cmap(:,:,4)=255;
cmap=permute(repmat(cmap,[1 3 1]),[2 1 3]);
%Repeat color for the 3 arrow vertices
set(q.Head,'ColorBinding','interpolated','ColorData', reshape(cmap(1:3,:,:), [], 4).'); 
% Repeat color for 2 tail vertices
set(q.Tail,'ColorBinding','interpolated','ColorData', reshape(cmap(1:2,:,:), [], 4).');
