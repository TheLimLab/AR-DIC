function out_frame=format_frame(frame,tile_size,scale)
% FORMAT_FRAME Prepares displacement frames for tiling in Tangram mapping software.
%
% Usage: out_frame=format_frame(frame,tile_size,scale)
%
% Input: frame is nxm displacement array. tile_size (scalar) sets pixel 
% dimensions of frame. A standard Tangram tile size is 256x256 pixels. 
% scale (scalar) sets scaling factor of frame. If needed, pads  the frame 
% with zeros to fill up to tile_size. Applies flipud to flip matrix 
% vertically.
%
% Returns: formatted displacement array.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% EXPORT_TO_MAP

temp_zero=zeros(tile_size,tile_size);%get zeros for padding
Z_temp=imresize(frame,scale); %scale
Z_temp=Z_temp.*scale;
Z_temp=flipud(Z_temp); %flip to match convention of other plots
[newZr,newZc]=size(Z_temp);
temp_zero(1:newZr,1:newZc)=Z_temp; %pad with zeros
out_frame=temp_zero;
end