function export_tree_centroids(leaves,obj,filename)
% EXPORT_TREE_CENTROIDS Exports centroids of all contracting regions to csv.
% Starting with the leaves, the tree is traversed from each leaf to the
% root.
%
% Usage: export_tree_centroids(leaves,obj,filename)
%
% Input: leaves, 1xn cell array containing leaves of tree with each leaf
% being a nodeobj. obj, vidobj filename, (string) output filename to save
% to.
%
% Output:csv file. Each node is exported in form of x, y, t triplets. Each
% row in file is path from leaf to root of tree. For example a tree with 3
% nodes: node1, node2, node3, with centroids (x1,y1), (x2,y2), (x3,y3)
% respectively a file may have the form: 
% x1,y1,t1,x2,y2,t2
% x3,y3,t1,x2,y2,t2,x3,y3,t3
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% REGION_TREE
% VIDOBJ

[yind,~]=size(obj.x_mat);
num_leaves=length(leaves);%initialize cells for storage in loop
leaf_path={};
lines={};
for n=1:num_leaves %for each leaf:
    temp_vect=[];
    leaf_path=trace_path(leaves{n},'time_sort');%get path from leaf to root
    for m=1:length(leaf_path)%for each node in path
       morph_ID=leaf_path{m}.ID(2);
       temp_cent=leaf_path{m}.data.morphology_props(morph_ID).Centroid; %get centroid and format it
       temp_cent(2)=yind-temp_cent(2);
       temp_vect(end+1)=temp_cent(1);
       temp_vect(end+1)=temp_cent(2);
       temp_vect(end+1)=leaf_path{m}.ID(1);%store centroid values and node ID    
    end
    lines{n}=temp_vect;
end
dlmwrite([filename '.csv'],lines{1})%write line to file
for c=2:length(lines)
dlmwrite([filename '.csv'],lines{c},'-append') %append additional lines to file
end

