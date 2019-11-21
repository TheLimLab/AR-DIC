function play_animation(Movie_structure,times,fps)
% PLAY_ANIMATION code is based almost line for line on Mathwork's example
% in 'Movie' Help document code for resizing animation frames and playing
% back with proper axis.
%
% Usage: play_animation(Movie_structure,times,fps)
%
% Inputs: Movie_structure is the Matlab movie object, times (scalar)
% specifies the number of times to play the movie. fps (scalar) specifies
% the frame rate.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln

F=Movie_structure;
% use 1st frame to get dimensions
[h, w, ~] = size(F(1).cdata);
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf,'Position', [150 150 w h]);
axis off
% Place frames at bottom left
movie(hf,F,times,fps,[0 0 0 0]);