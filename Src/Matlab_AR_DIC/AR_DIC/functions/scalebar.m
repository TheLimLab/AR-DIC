function scalebar(xstart,ystart,W,H,color)
% SCALEBAR Draws a scale bar rectangle over current image.
%
% Useage: scalebar(xstart,ystart,W,H,color)
%
% Input: xstart and ystart are scalars that define the horizontal and
% vertical coordinates for rectangle origin. W, H are scalars which define
% the width and height of the rectangle. color defines the scale bar color
% using any of the default Matlab facecolor options.
%     
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln

rectangle('pos',[xstart ystart W H],'facecolor',color,'edgecolor','none')