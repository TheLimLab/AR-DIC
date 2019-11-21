classdef nodeobj < handle
% NODEOBJ stores data and points to children and parent nodes. 
% Suitable for building data trees. Tracks both parent and child node to
% facilitate top-down building from object region_tree. If constructed with
% no inputs then ID is set to 0. nodeobj is called automatically from
% region_tree to build a data tree.
%
% Usage: obj=nodeobj(parent, data, ID)
%
% Inputs: (variable number of inputs): parent (nodeobj), data (frameobj),
% ID ([frame_number, morphology_index]).
%
% NODEOBJ Properties:
%   data - Stores data here
%   children - Cell array pointing to children handles, if node has no children then cell array is empty
%   ID - Node identification, of form [frame_number, morphology_index]
%   parent - Points to parent node
%   parent_ID - Identification of parent
%
% NODEOBJ Methods:
%   nodeobj - Object constructor. varargin order: parent, data, ID
%   trace_path - Traces path from the input node obj to root
%   plot_region - Plots binary mask region associated with node
%
% Object hierarchy:
% 1. vidobj
%   i. frameobj
%  ii. region_tree
%       a. NODEOBJ (has access to frameobj data)
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% REGION_TREE
% FRAMEOBJ
% VIDOBJ
    
    properties
        data %Stores data here
        children %Cell array pointing to children handles, if node has no children then cell array is empty
        ID %Node identification, of form [frame_number, morphology_index]
        parent %Points to parent node
        parent_ID %Identification of parent
    end
    methods
        function obj=nodeobj(varargin)
        %NODEOBJ Constructor.
        %nodeobj is intended to be called automatically from region_tree
        %which constructs a data tree and creates nodeobjs on the fly.
        %
        %Usage: obj=nodeobj(varargin)
        %
        %Input: varargin order: parent, data, ID.
        %
        %Returns: a new nodeobj
            
            
            if nargin==0 %if constructed with no input, create empty object
                obj.ID=0; %root node
                obj.children={};
            elseif nargin==1 %if only parent is given
                obj.parent=varargin{1};
                obj.children={};
            else %if parent and data are given
                obj.parent=varargin{1};
                obj.data=varargin{2};
                obj.children={};
            end
        end% nodeobj constructor
        function out_array=trace_path(obj,sorttype)
        %TRACE_PATH Traces path from the input node obj to root. 
        %
        %Usage: out_array=trace_path(obj,sorttype)
        %
        %Input: sorttype (string) specifies which sort method to use.
        %Current option: 'time_sort'.
        %
        %Returns: out_array, a 1xn cell array containing nodeobjs in path
        %order from obj to root of data tree.
            
            out_array{1}=obj;%initialize output array starting with leaf
            nodetype=class(obj);%get class of obj.
            parent_node=obj.parent;%get parent of obj.
            while isa(parent_node,nodetype) %continue adding nodes until root is reached, at which point parent node is a vidobj instead of nodeobj
                out_array{end+1}=parent_node; %store node sequence
                parent_node=parent_node.parent; %update parent
            end
            switch sorttype
                case 'time_sort'
                    out_array=sort_path_time(out_array);%sort based on time in ascending order
            end
        end%trace_path
        function plot_region(obj,mag)
        %PLOT_REGION Plots binary mask region associated with node.
        %Contracting regions are plotted in white on black background.
        %
        %Usage: plot_region(obj,mag)
        %
        %Input: mag (scalar) sets initial image magnification.
            
            [ydim,xdim]=size(obj.data.xdisp);
            L=zeros(ydim,xdim);%initialize background
            data_ind=obj.ID(2);
            CC=obj.data.CC.PixelIdxList{data_ind};%get contracting objects
            L(CC)=1;%set contracting objects to 1
            imshow(L,'InitialMagnification',mag)%display with magnification defined by mag
            axis xy %match convention
            axis([0 xdim 0 ydim])%set axis limits
        end
    end%methods
end %nodeobj

%utility function:
function time_nodes=sort_path_time(path_cell)
%SORT_PATH_TIME Sorts cell array of nodeobjs based on frame number (obj.ID(1)).
%Function is called from trace_path to sort node path in ascending order
%based on time.
%
%Usage: time_nodes=sort_path_time(path_cell)
%
%Input: path_cell (1xn cell array) containing nodeobjs to be sorted.
%
%Returns: sorted 1xn cell array of nodeobjs.
%
%See also:
%NODEOBJ/TRACE_PATH

num_nodes=length(path_cell);
times=zeros(1,num_nodes);
for n=1:num_nodes%get times
    times(n)=path_cell{n}.ID(1);%get time from first index of id
end
[~,index]=sort(times,'ascend'); %sort time to get index
time_nodes=path_cell(index);%sort nodes in descending order based on index
end%sort_path_time