classdef frameobj < handle
% FRAMEOBJ Is a handle class for managing and analyzing DIC mechanics data.
% In general, frameobjs are intended to be created automatically through
% the vidobj. Each frame of the vidobj is organized into a frameobj.
%
% Usage: obj=frameobj(varargin)
% 
% Inputs: filename (string), scale (scalar), timestep (scalar), xprev (nxm
% array of x-displacements) , yprev (nxm array of y-displacements,
% refmethod (optional, currently unused input to set reference scoring
% method).
%
% FRAMEOBJ Properties:
%  binary_mask - Binary mask from threshold
%  CC - Connected components in binary frame from bwconncomp
%  contraction_volume - Displacement magnitude times thresholded area
%  disp_mat - Displacement magnitude matrix
%  disp_max_index - Index of maximum displacement
%  disp_min_index - Index of minimum displacement
%  file - Input filename and paths
%  max_disp - Maximum displacement in frame
%  max_location_x - Maximum displacement x location
%  max_location_y - Maximum displacement y location
%  max_vel - Maximum velocity in frame
%  min_disp - Minimum displacement in frame
%  min_location_x - Minimum displacement x location
%  min_location_y - Minimum displacement y location
%  morphology_props - Morphology measurement structure
%       note: Centroid, BoundingBox, and ConvexHull are stored unscaled.
%  percent_area - Percent of frame above displacement threshold
%  refscore - refscore is currently unused in frameobj
%  strain_tensor - Structure containing strain tensor
%  u_mat - x direction velocity matrix
%  v_mat - y direction velocity matrix
%  vel_mat - Velocity magnitude matrix
%  xdisp - x direction displacement matrix
%  ydisp - y direction displacement matrix
%
% FRAMEOBJ Methods:
%  calc_contraction_volume - Calculates contraction volume for frame and sets in object properties
%  calc_morphology - Process morphology data from binary masks and sets morphology properties of frameobj 
%  calc_principal_strain - Calculate principal strains from raw strain data
%  calc_tensor - Calculates strain tensor and sets in object properties
%  find_mask - Calculate binary mask based on threshold
%  frameobj - Object constructor
%  get_max_morphology - Retrieves maximum value from morphology data
% 
% Object hierarchy:
% 1. vidobj
%   i. FRAMEOBJ
%  ii. region_tree
%       a. nodobj (has access to frameobj data)
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
% 
% See also:
% VIDOBJ
% NODEOBJ
% REGION_TREE

properties
    binary_mask %Binary mask from threshold
    CC %Connected components in binary frame from bwconncomp
    contraction_volume %Displacement magnitude times thresholded area
    disp_mat %displacement magnitude matrix
    disp_max_index %Index of maximum displacement
    disp_min_index %Index of minimum displacement
    file %Input filename and paths
    max_disp %Maximum displacement in frame
    max_location_x %Maximum displacement x location
    max_location_y %Maximum displacement y location
    max_vel %Maximum velocity in frame
    min_disp %Minimum displacement in frame
    min_location_x %Minimum displacement x location
    min_location_y %Minimum displacement y location
    morphology_props %Morphology measurement structure
%       note: Centroid, BoundingBox, and ConvexHull are stored unscaled.
    percent_area %Percent of frame above displacement threshold
    refscore %refscore is currently unused in frameobj
    strain_tensor %Structure containing strain tensor
    u_mat %x direction velocity matrix
    v_mat %y direction velocity matrix
    vel_mat %Velocity magnitude matrix
    xdisp %x direction displacement matrix
    ydisp %y direction displacement matrix
end
methods
    function obj=frameobj(varargin) 
    %FRAMEOBJ constructor
    %
    %Usage: obj=frameobj(varargin)
    %
    %Input order: filename, scale, timestep, xprev(optional),
    %yprev(optional), method(optional). frameobj can be created with either
    %0 inputs, the first 3 inputs, or all inputs.
    %
    %Returns: a new frameobj.
        
            if nargin==0 %if no input create empty object
                obj.file={};
            else
                filename=varargin{1};%record data file name
                scale=varargin{2};
                timestep=varargin{3};
                obj.file=filename;
                if nargin==3
                %Read data from DIC text file
                    [x_mat,y_mat,obj.xdisp,obj.ydisp,obj.disp_mat,obj.u_mat,...
                        obj.v_mat,obj.vel_mat]=readimagejPIV(filename,scale,timestep);
                else
                    xprev=varargin{4};
                    yprev=varargin{5};
                    refmethod=varargin{6};
                    %Read data from DIC text file
                    [x_mat,y_mat,obj.xdisp,obj.ydisp,obj.disp_mat,obj.u_mat,...
                        obj.v_mat, obj.vel_mat]=readimagejPIV(filename,scale,timestep,xprev,yprev); %Format data
%                     obj.refscore=calcrefscore(obj.disp_mat,refmethod);
                end
                %find maximum and minimum values, store in object
                [disp_max,maxindex]=max(obj.disp_mat(:));%find maximum displacement magnitude and location
                [I_row, I_col]=ind2sub(size(obj.disp_mat),maxindex); %get row, column indices from maxindex
                x_maxlocation=x_mat(I_row,I_col);
                y_maxlocation=y_mat(I_row,I_col);
                obj.disp_max_index=[I_row,I_col];
                obj.max_disp=disp_max; %store min disp in object
                obj.max_location_x=x_maxlocation; %store x, y in object
                obj.max_location_y=y_maxlocation;
                [disp_min,minindex]=min(obj.disp_mat(:));%find maximum displacement magnitude and location
                [I_row_min, I_col_min]=ind2sub(size(obj.disp_mat),minindex);
                x_minlocation=x_mat(I_row_min,I_col_min);
                y_minlocation=y_mat(I_row_min,I_col_min);
                obj.disp_min_index=[I_row_min,I_col_min];
                obj.min_disp=disp_min; %store min disp in object
                obj.min_location_x=x_minlocation; %store x,y in object
                obj.min_location_y=y_minlocation;
                obj.max_vel=max(obj.vel_mat(:));%max velocity magnitude
            end
        end %end class constructor
        function find_mask(obj,threshold)
        %FIND_MASK Calculates binary mask based on threshold.
        %Values above the threshold are assigned a 1 and everything else is
        %filled with 0.
        %
        %Usage: find_mask(obj,threshold)
        %
        %Input: threshold (scalar)
        
            one_frame=obj.disp_mat;
            obj.binary_mask=one_frame>threshold;
        end
        function calc_contraction_volume(obj)
        %CALC_CONTRACTION_VOLUME Calculates contraction volume for frame and sets in object properties
        %
        %Usage: calc_contraction_volume(obj)
        
            thresholded=obj.disp_mat.*obj.binary_mask;
            obj.contraction_volume=sum(thresholded(:));
        end
        function calc_tensor(obj,h,type)
        %CALC_TENSOR Calculates strain tensor and sets in object properties.
        %
        %Usage: calc_tensor(obj,h,type)
        %
        %Inputs: h (scalar) is spacing between sample points. type (string)
        %determines calculation method: 'cauchy' or 'green_lagrangian'.
        
            ux=obj.xdisp;
            uy=obj.ydisp;
            [Uxx,Uxy]=gradient(ux,h);
            [Uyx,Uyy]=gradient(uy,h);
            switch type
                case 'cauchy'
                    E_xx=Uxx;
                    E_yy=Uyy;
                    E_xy=(1/2)*(Uxy+Uyx);
                case 'green_lagrangian'
                    E_xx=Uxx+0.5*((Uxx.^2)+(Uyx.^2));
                    E_yy=Uyy+0.5*((Uxy.^2)+(Uyy.^2));
                    E_xy=0.5*(Uxy+Uyx)+0.5*((Uxx.*Uxy)+(Uyx.*Uyy));
            end  
            obj.strain_tensor.Exx=E_xx;
            obj.strain_tensor.Eyy=E_yy;
            obj.strain_tensor.Exy=E_xy;            
        end%calc_tensor
        function out=calc_principal_strain(obj)
        %CALC_PRINCIPAL_STRAIN Calculates principal strains from raw strain data. 
        %
        %Usage: out=calc_principal_strain(obj)
        %
        %Returns: out, structure with fields s1, s2, theta, and
        %max_eng_shear containing matrices of maximum principal strain,
        %minimum principal strain, transform angle, and maximum engineering
        %strain respectively.
        
            E_xx=obj.strain_tensor.Exx;%get raw strain matrices
            E_yy=obj.strain_tensor.Eyy;
            gammaxy=2*(obj.strain_tensor.Exy);
            etam=(E_xx+E_yy)/2;
            etadiff=(E_xx-E_yy)/2;
            R=sqrt(etadiff.^2+gammaxy.^2);
            p_strain1=etam+R; %principal strain component1
            p_strain2=etam-R; %principal strain component2
            theta_p_strain=(1/2).*atand(gammaxy./(E_xx-E_yy)); %angle of rotation [degrees] for principal strain
            out.s1=p_strain1;%set principal strain values
            out.s2=p_strain2;
            out.theta=theta_p_strain;
            out.max_eng_shear=R;
        end
        function calc_morphology(obj,scale,type)
        %CALC_MORPHOLOGY Process morphology data from binary masks and sets morphology properties of frameobj
        %vidobj/measure_morphology calls calc_morphology for each frame in
        %vidobj
        %
        %Usage: calc_morphology(obj,scale,type)
        %
        %Inputs: scale (scalar) is the pixel scale (e.g. px/um). type
        %(string) either 'box' or 'convex_hull'. Specifies if BoundingBox
        %or ConvexHull should be calculated by regionprops.
        %
        %See also:
        %VIDOBJ/MEASURE_MORPHOLOGY
        
        switch type
            case 'box'
            props={'Area', 'Perimeter', 'MajorAxisLength',...
                'MinorAxisLength','Eccentricity','Orientation',...
                'EquivDiameter','Centroid','BoundingBox'};%area properties to measure
            case 'convex_hull'
            props={'Area', 'Perimeter', 'MajorAxisLength',...
                'MinorAxisLength','Eccentricity','Orientation',...
                'EquivDiameter','Centroid','ConvexHull'};%area properties to measure
        end
            C_region=bwconncomp(obj.binary_mask,8); %Find connected regions
            morph_props=regionprops(C_region,props); %compute region properties
            obj.morphology_props=scale_morphology(morph_props,scale); %use scale on morphology
            obj.percent_area=100*(sum(obj.binary_mask(:))/length(obj.binary_mask(:)));
            obj.CC=C_region; %store connected components from bwconncomp
        end
        function max_prop=get_max_morphology(obj,morph_prop)
        %GET_MAX_MORPHOLOGY Retrieves maximum value of 'morph_prop' from morphology data.
        %
        %Usage: max_prop=get_max_morphology(obj,morph_prop)
        %
        %Input: morph_prop (string) specifies property to measure. Any of
        %the morphology structure parameters: 'Area', 'Perimeter',
        %'MajorAxisLength','MinorAxisLength','Eccentricity','Orientation',
        %'EquivDiameter','Centroid'.
        %
        %Returns: max_prop (scalar), value of maximum morphology value.
        
            frame1=obj.morphology_props;
            temp_max=0; %initialize max to 0
            for n=1:length(frame1) %collect maximums and find greatest
                temp_area=frame1(n).(morph_prop);
                if temp_area>temp_max
                    temp_max=temp_area;
                end
            end
            max_prop=temp_max;
        end%get_max_morphology
    end %methods end
end %classdef end
function scaled_st=scale_morphology(input_st,scale)
%SCALE_MORPHOLOGY Called automatically from calc_morphology. 
%Scales morphology measurements: Area, major axis length, minor axis
%length, equivalent diameter, and perimeter by scale.
%
%Usage: scaled_st=scale_morphology(input_st,scale)
%
%Inputs: input_st is 1xn output of regionprops function. scale (scalar)
%defines scale to correct morphology measurements.
%
%Returns: scaled_st, 1xn array of scaled morphology values
%
%See also:
%FRAMEOBJ/CALC_MORPHOLOGY

scaled_st=input_st;
    for n=1:length(input_st)
        scaled_st(n).Area=input_st(n).Area*scale^2;
        scaled_st(n).MajorAxisLength=input_st(n).MajorAxisLength*scale;
        scaled_st(n).MinorAxisLength=input_st(n).MinorAxisLength*scale;
        scaled_st(n).EquivDiameter=input_st(n).EquivDiameter*scale;
        scaled_st(n).Perimeter=input_st(n).Perimeter*scale;
    end
end