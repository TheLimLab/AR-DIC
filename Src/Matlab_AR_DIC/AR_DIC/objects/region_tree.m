classdef region_tree < handle
% REGION_TREE Constructs data tree by organizing segmented areas from vidobj. 
% Expects to use binary masks obtained from vidobj method mask_threshold.
% Establishes parent/child relationship between all areas.
%
% Usage: obj=region_tree(in_vidobj,type)
%
% Input: in_vidobj is a vidobj, type (string) either 'box' or
% 'convex_hull'. Specifies method to use to check if child bounded by
% parent region. 'box' may be significantly faster but may not deal with
% complex shapes well. 'convex_hull' is more accurate but may be
% significantly slower.
%
% REGION_TREE Properties:
%  leaf_index - (Dependent property) Stores index to every node in the region_tree
%  node_holder - Stores tree structure with linked nodes
%  root - For bookkeeping, root of tree is the video (a vidobj)
%
% REGION_TREE Methods:
%  ASC_map - Visualize region tree in either 2D or 3D
%  check_lineage - Finds parent of input node, assigns parent-child data, returns level of child node
%  get_node - Retrieves node(s) from region_tree at 2xn array node_index
%  region_tree - Object constructor
%
% Object hierarchy:
% 1. vidobj
%   i. frameobj
%  ii. REGION_TREE
%       a. nodobj (has access to frameobj data)
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% NODEOBJ
% FRAMEOBJ
% VIDOBJ
% VIDOBJ/MASK_THRESHOLD

    properties
        node_holder %Stores tree structure with linked nodes
        root %For bookkeeping, root of tree is the video (a vidobj)
    end
    properties (Dependent = true)
        %Stores index to every node in the region_tree. Array: [level_index,cell_index]
        leaf_index 
    end
    methods
        function obj=region_tree(in_vidobj,type) 
        %REGION_TREE Constructor. Builds data tree from areas contained in vidobj.
        %Establishes parent/child relationship between all areas. Expects
        %to use binary masks obtained from vidobj method mask_threshold. A
        %region_tree is constructed automatically when construct_tree is
        %called on a vidobj.
        %
        %Usage: obj=region_tree(in_vidobj,type)
        %
        %Input: in_vidobj is a vidobj, type (string) either 'box' or
        %'convex_hull' specifies method to use to check if child bounded by
        %parent region. 'box' may be significantly faster but may not deal
        %with complex shapes well. 'convex_hull' is more accurate but may
        %be significantly slower. Future releases will add more options.
        %
        %Returns: a new region_tree object
        
            obj.root=in_vidobj; %set input vidobj as root 
            sorted=sort_all_areas(in_vidobj); %sort all areas by descending size
            [~,num_areas]=size(sorted);
            obj.root=in_vidobj;%make parent node
            first=nodeobj(obj.root);%make first node from largest area
            first.parent_ID=0;
            first.ID=sorted(2:3,1); %set ID with frame and morphology index
            first.data=in_vidobj.frame_holder{1,first.ID(1,1)}; %point to data
            obj.node_holder{1}={first};
            for n=2:num_areas %iterate through all areas to build tree
                temp=nodeobj();%new node
                temp.ID=sorted(2:3,n); %set node ID
                temp.data=in_vidobj.frame_holder{1,temp.ID(1,1)}; %point to data
                [level]=check_lineage(obj,temp,type);%finds parent of node, assigns parent, gives child level
                if level>length(obj.node_holder)
                    obj.node_holder{level}={temp};%if node is in higher level, add node
                else
                    obj.node_holder{level}{end+1}=temp;%otherwise add node to current level
                end
            end
        end %region_tree constructor
        function out=get.leaf_index(obj)
        %GET.LEAF_INDEX Finds all leaves of data tree and return array of indices to leaves. 
        %
        %Usage: out=get.leaf_index(obj)
        %
        %Returns: [level_index;branch_index] which provides the full path to each leaf.
        
            level_index=[];
            branch_index=[];
            for nL=1:length(obj.node_holder)%iterate through each node
                level=obj.node_holder{1,nL};%get level
                for n=1:length(level)%find leaves
                    if isempty(level{n}.children)%if node has no children then it is a leaf
                        level_index(end+1)=nL;
                        branch_index(end+1)=n;
                    end
                end
            end
            out=[level_index;branch_index];%return indices to leaves
        end %leaf_index
        function out=get_node(obj,node_index)
        %GET_NODE Retrieves node(s) from region_tree at 2xn array node_index.
        %
        %Usage: out=get_node(obj,node_index)
        %
        %Input: node_index has the form [level_index;cell_index].
        %
        %Returns: 1xn cell array containing specified nodeobjs.
        
            [~,num_leaves]=size(node_index);
            out=cell(1,num_leaves);
            for n=1:num_leaves
                level_ind=node_index(1,n);
                cell_ind=node_index(2,n);
                out{n}=obj.node_holder{level_ind}{cell_ind};
            end
        end %leaf_array
        function ASC_map(obj,view_type,pause_time)
        %ASC_MAP Visualizes accumulative spatial contraction map in either 2D or 3D. 
        %The 3D option produces the Accumulative Spatial Contraction (ASC)
        %map.
        %
        %Usage: ASC_map(obj,view_type,pause_time)
        %
        %Inputs: view_type (string) can be '2D' or '3D' for 2D or 3D
        %visualization respectively. pause_time (scalar) and is for
        %animation purposes only and can be set to zero.
        
            AZ=-37.5;%default: -37.5
            EL=30; %default: 30
            [ydim,xdim]=size(obj.root.x_mat);
            zmax=length(obj.root.frame_holder);
            num_level=length(obj.node_holder);
            L=zeros(ydim,xdim);
            figure
            for n=1:num_level%iterate through levels, plot regions sequentially
                for n2=1:length(obj.node_holder{1,n})
                    data_ind=obj.node_holder{1,n}{1,n2}.ID(2);%get index to each region
                    CC=obj.node_holder{1,n}{1,n2}.data.CC.PixelIdxList{data_ind};
                    switch view_type
                        case '2D'%2D view of ASC map
                            L(CC)=n;%new label for each region
                            H=label2rgb(L);%convert label to color
                            imshow(H,'InitialMagnification',600)
                            axis xy
                            axis([0 xdim 0 ydim])
                        case '3D' %Accumulative Spatial Contraction map
                            L(CC)=n;%set value to current level
                            surf(obj.root.x_mat,obj.root.y_mat,L)
                            view(AZ,EL)
                            axis ij
                            colorbar
                    end
                    pause(pause_time)%pause time only used for visualization purposes
                end
                hold on
            end
        end%ASC_map
        function [level]=check_lineage(obj,test_node,type)
        %CHECK_LINEAGE Given a region_tree object, finds parent of nodeobj test_node.
        %Also assigns parent-child data and returns level of child node. 
        %
        %Usage: [level]=check_lineage(obj,test_node,type)
        %
        %Input: test_node, the nodeobj to find parents for, type (string)
        %either 'box' or 'convex_hull' which specifies which method to use
        %to check if the child is bounded by the parent region. 'box' may
        %be significantly faster but may not deal with complex shapes well.
        %'convex_hull' is more accurate but may be significantly slower.
        %Future releases will add more options.
        %
        %Returns: level (scalar) of child node.
        
            morph_id=test_node.ID(2);%get node morphology ID
            cent=test_node.data.morphology_props(morph_id).Centroid;%get centroid from test node
            levels=length(obj.node_holder);
            for lev=levels:-1:1 %countdown loop
                for n=1:length(obj.node_holder{lev})
                    test_parent=obj.node_holder{lev}{n};%set test parent
                    pm_id=test_parent.ID(2);%get node morphology ID
                    switch type
                        case 'box'
                        box=test_parent.data.morphology_props(pm_id).BoundingBox; %get parent bounding box
                        case 'convex_hull'
                        box=test_parent.data.morphology_props(pm_id).ConvexHull;
                    end
                    isparent=in_region(box,cent,type); %check if centroid is in bounding box    
                    if isparent %if centroid is in bounding box, test_parent is parent
                        test_node.parent=test_parent; %set parent
                        test_node.parent_ID=test_parent.ID;%set parent id
                        test_parent.children{end+1}=test_node;%set children of parent
                        break
                    end
                end
                if isparent
                    break
                end
            end
            if ~isparent %if parent is not in tree, then it is an independent region, root is parent
                test_node.parent=obj.root;%set parent as root
                test_node.parent_ID=0;%set root id
                level=1; %since parent is root, child is in first level
            else
                level=lev+1; %increment level from parent level to store child node
            end
        end%check_lineage
    end %methods
end%region_tree
function binary_out=in_region(bound_region,cent,type)
%IN_REGION Checks if centroid, cent, is bounded by region bound_region.
%Output true if centroid is in bounding box, otherwise output is false.
%
%Usage: binary_out=in_region(bound_region,cent,type)
%
%Input: bound_region is either a 1x4 array as returned from regionprops
%BoundingBox property or the convex hull returned by the regionprops
%ConvexHull property. BoundingBox format: [upperleft x-coordinate,
%upperleft y-coordinate, x-width, y-width]. cent (1x2 array) containing the
%[x,y] coordinates of the centroid to test. type (string) either 'box' or
%'convex_hull'. Specifies method to use to check if child bounded by parent
%region. 'box' may be significantly faster but may not deal with complex
%shapes well. 'convex_hull' is more accurate but may be significantly
%slower. Future releases will add more options.
%
%Returns: logical 'true' if centroid is bounded by region otherwise returns
%'false'.

x=cent(1);
y=cent(2);
switch type
    case 'box' %use box method to check if centroid (x,y) is in bounding box
        ulx=bound_region(1);
        uly=bound_region(2);
        x_width=bound_region(3);
        y_width=bound_region(4);
        inx_region=(x>ulx)&&(x<(ulx+x_width));
        iny_region=(y>uly)&&(y<(uly+y_width));
        if inx_region&&iny_region
            binary_out=true;
        else
            binary_out=false;
        end
    case 'convex_hull'%use convex hull to check if centroid (x,y) is inside convex hull
        polyg_x=bound_region(:,1);
        polyg_y=bound_region(:,2);
        binary_out=inpolygon(x,y,polyg_x,polyg_y);
end
end%in_region