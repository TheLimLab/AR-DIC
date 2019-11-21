function plot_node_array(leaf_array,obj)
% PLOT_NODE_ARRAY Plots centroids of nodes contained in leaf_array.
%
% Usage: plot_node_array(leaf_array,obj)
%
% Inputs: leaf_array is a 1xn cell array of nodeobj containing the leaf
% nodes. obj is the associated vidobj.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% NODEOBJ

[yind,~]=size(obj.x_mat);%get dimension
figure
hold on
for n=1:length(leaf_array)%for each node in array
    data_ind=leaf_array{1,n}.ID(2);%get data location
    cent=leaf_array{1,n}.data.morphology_props(data_ind).('Centroid'); %get centroid
    cent(2)=yind-cent(2); %format centroid
    cent=cent.*obj.vect_space.*obj.scale;%use scale factors for plotting
    scatter(cent(1),cent(2),60,'filled','k') %plot centroids
end
hold off
title('Smallest independent regions')
end