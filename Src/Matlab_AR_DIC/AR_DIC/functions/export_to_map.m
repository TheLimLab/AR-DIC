function [RGBmat,Alpha]=export_to_map(vidobj_holder,start_stop,num_col,high_alpha,low_alpha,tile_size,opt)
% EXPORT_TO_MAP Transforms and tiles displacement frames. 
% The output format is suitable for viewing in open source map rendering
% software Tangram. The displacement is encoded in the surface normals of
% an RGB plot and height data is translated to the alpha channel.
%
% Usage: [RGBmat,Alpha]=export_to_map(vidobj_holder,start_stop,num_col,high_alpha,low_alpha,tile_size,opt)
%
% Input: vidobj_holder is a 1xn cell array which holds vidobjs to plot. The
% function loops through each vidobj and tiles the specified frames from
% each. start_stop is a 1xn cell array which defines the range of frames to
% tile. Each element of the cell array is of the format [start_index,
% stop_index]. num_col is a scalar which defines the number of columns to
% divide the specified frames into. num_col must be defined so that the
% number of rows times the number of columns equals the total number of
% frames. high_alpha is a scalar which defines the upper limit of the
% height map. low_alpha is a scalar which defines the lower range of the
% height map. Using the encoding scheme in Tangram sea level has an alpha
% value of 0.93. Mt. Diablo in CA, USA has an alpha height of 0.82.
% tile_size is a scalar which defines the size of the tile block. A
% standard Tangram tile is 256x256 pixels. The displacement matrix is
% scaled to fit the tile size. If the displacement frame is less than the
% tile_size then the block is padded with zeros. opt is a string, either
% 'maxmin' to plot only the maximum and minimum displacement frames or
% 'all' to plot the frames specified by start_stop. Returns: RGBmat an
% nxmx3 matrix containing the tiled displacements with all surface normal
% values scaled from 0 to 255. Alpha is an nxmx3 matrix containing the
% height data scaled to the alpha channel.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% FORMAT_FRAME

[tempr,tempc]=size(vidobj_holder{1}.frame_holder{1}.disp_mat); %get size of displacement matrix
maxdim=max([tempr,tempc]);
scale=floor(tile_size/maxdim); %get scale to transform to tile size
num_videos=length(vidobj_holder);
switch opt
    case 'maxmin' %get max and min stretch frames
        Z=[];
        for nt=1:num_videos%loop through all videos
            vid=vidobj_holder{nt};
            temp_max=vid.frame_holder{1,vid.max_disp_frame}.disp_mat; %maximum displacement frame
            temp_min=vid.frame_holder{1,vid.min_disp_frame}.disp_mat; %minimum displacement frame
            %Pad frames with zeros, scale, format, and concatenate for
            %later processing
            format_max=format_frame(temp_max,tile_size,scale);
            format_min=format_frame(temp_min,tile_size,scale);
            Z=cat(2,Z,[format_min;format_max]); %organize frames
        end
        
    case 'all' %get all frames
        Z=[];
        for vids=1:length(vidobj_holder)
            counter=start_stop{vids}; %get frame index and format
            start_frame=counter(1);
            stop_frame=counter(2);
            n_count=start_frame:stop_frame;
            num_row=length(n_count)/num_col;
            temp_zero=cell(num_row,num_col);
            n=1;
            for row=1:num_row
                for col=1:num_col
                matmat=vidobj_holder{vids}.frame_holder{1,n_count(n)}.disp_mat;%get displacement frame
                temp_zero{row,col}=format_frame(matmat,tile_size,scale);%format and store frame
                n=n+1;
                end
            end
            hold_mat=cell2mat(temp_zero);%convert cell array to matrix for later processing
            Z=cat(2,Z,hold_mat); %organize frames
        end
end
minZ=min(Z(:)); %get max and min to scale alpha elevation
maxZ=max(Z(:));
[R,G,B]=surfnorm(Z); %get surface normals
Alpha=(Z-minZ)*(high_alpha-low_alpha)/(maxZ-minZ)+low_alpha;% scale alpha elevation
RGBmat=reshape([R G B], size(R,1), size(R,2),3); %place RGB data in matrix format
RGBmat=mat2gray(RGBmat); %scale to 0 to 255 range
