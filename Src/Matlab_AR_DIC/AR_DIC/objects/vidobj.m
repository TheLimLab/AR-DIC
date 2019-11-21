classdef vidobj < handle
% VIDOBJ is a handle class that processes and stores data from digital image correlation files.
% vidobj uses stringtext and folderpath to open DIC output files
% sequentially. vidobj creates a frameobj for each DIC file.
%
% Usage: obj=vidobj(stringtext, folderpath, scale, timestep). 
% 
% Inputs: stringtext (string) specifies the text to search for to input DIC
% files. For example, for files from ImageJ PIV plugin set stringtext to
% *PIV3*. vidobj sorts thse filenames in numerical order for processing.
% folderpath (string) specifies the directory containing DIC files. scale
% (scalar) specifies the pixel scale (e.g. px/um). timestep (scalar)
% specifies the time between image frames.
%
% VIDOBJ Properties:
% x_mat - x location matrix
% y_mat - y location matrix
% xdim - x dimension of frame
% ydim - y dimension of frame
% num_frame - Number of frames
% scale - Scale of frame (e.g. px/um)
% vect_space - Spacing of DIC sample vector
% timestep - Time between frames
% frame_holder - Cell array for holding frameobjs
% max_disp - Maximum displacement from all frames
% max_disp_frame - Frame number with maximum displacement
% min_disp - Minimum displacement from all frames
% min_disp_frame - Frame number with minimum displacement
% max_velocity - Maximum velocity from all frames
% max_velocity_frame - Frame with maximum velocity
% tree - Tree data structure
%
% VIDOBJ Methods:
% vidobj - Object constructor
% Plotting methods
% plot_contour - Creates either displacement or velocity contour plot
% plot_quiver - Creates either displacement or velocity quiver plot
% plot_threshold_contour - Creates contour plot of displacement values above specified threshold
% heatmap_frequency - Creates heatmap of contraction frequency
% heatmap_magnitude - Creates heatmap of contraction magnitude
% surfplot - Visualize displacement in 3D plot using z axis to encode magnitude
% plot_contraction_volume - Plots the contraction volume vs. time
% plot_binary_mask - Animation formed by plotting each thresholded binary mask of vidobj
% xy_quiver - Animates quiver plot with displacement in x direction plotted with blue arrows and displacement in y direction with red arrows
% plot_centroid_areathreshold - Plots centroids for areas greater than specified area threshold
% plot_max_disp_each_frame - Plots the maximum displacement in each frame vs. time
% plot_disp_location - Plots displacement for each frame at location specified by index
% 
% Calculation methods
% percent_area_disp_thresh - Thresholds displacement to obtain percent contracting area above threshold
% mask_threshold - Obtains binary mask for contracting areas with displacements greater than threshold
% contraction_volume - Calculates the contraction volume for each frame
% measure_morphology - Calculates morphology properties by passing each frame to calc_morphology
% calc_strain - Calculates strain tensor for all video frames
% transform_strain - Transforms the raw strain data using the specified method
% strain_to_cell - Collects the strain tensor from each frameobj and returns in a cell array
% sort_morphology - Sorts the morphology measurements from largest to smallest
% sort_all_areas - Sorts all areas in descending order
% construct_tree - Passes vidobj to region_tree which constructs a data tree from the areas
% get_disp - Obtains displacement at index for every frame
% get_contraction_volume - Get contraction volume from all frames 
% get_max_velocity - Get maximum velocity from each frame
%
% Object hierarchy:
% 1. VIDOBJ
%   i. frameobj
%  ii. region_tree
%       a. nodobj (has access to frameobj data)
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% frameobj
% nodeobj
% region_tree

    properties
        x_mat %x location matrix
        y_mat %y location matrix
        xdim %x dimension of frame
        ydim %y dimension of frame
        num_frame %Number of frames
        scale %Scale of frame (e.g. px/um)
        vect_space %Spacing of DIC sample vector
        timestep %Time between frames
        frame_holder %Cell array for holding frameobjs
        max_disp %Maximum displacement from all frames
        max_disp_frame %Frame number with maximum displacement
        min_disp %Minimum displacement from all frames
        min_disp_frame %Frame number with minimum displacement
        max_velocity %Maximum velocity from all frames
        max_velocity_frame %Frame with maximum velocity
        tree %Tree data structure 
    end
    
    methods
        function obj=vidobj(varargin) 
        %VIDOBJ constructor
        %
        %Usage: obj=vidobj(varargin)
        %
        %Input: stringtext (string) specifies the text string to search for
        %in DIC filenames. For example, from ImageJ PIV plugin set
        %stringtext to '*PIV3*'. folderpath (string) specifies the
        %directory containing DIC files. scale (scalar) specifies the pixel
        %scale (e.g. px/um). timestep (scalar) specifies the time between
        %image frames.
        %
        %Returns: a new vidobj.
        
            if nargin==0%if no inputs create empty object
                obj.num_frame=0;
            else
                stringtext=varargin{1};
                folderpath=varargin{2};
                scale=varargin{3};
                timestep=varargin{4};
                directory=pwd; %get current directory
                cd(folderpath) %change directory to get input files
                PIVtext=dir(stringtext);%get DIC filenames
                PIVtext_names={PIVtext.name};
                [~,idx]=sort_nat(PIVtext_names);%sort filenames numerically
                PIVtext=PIVtext(idx);
                frames=length(PIVtext);
                %get dimension data from first frame
                [obj.x_mat,obj.y_mat,~,~,~,~,~,~]=readimagejPIV(PIVtext(1).name, scale, timestep);
                obj.xdim=obj.x_mat(1,end);
                obj.ydim=obj.y_mat(1,1);
                obj.vect_space=(obj.x_mat(1,2)-obj.x_mat(1,1)); %vector spacing after scaling by px/um ratio
                frame_hold=cell(1,frames);
                if nargin<5
                    for num=1:frames%format each frame into frameobj
                        filename=PIVtext(num).name;
                        frame_hold{num}=frameobj(filename, scale, timestep); %Format data
                    end
                else
                    %in future revisions measure ref score here 
                end
                %find and store maximum values
                val_max=zeros(1,frames); %preallocate storage
                val_min=val_max;%preallocate storage
                val_vel=val_max;%preallocate storage
                for num2=1:frames %maximum displacement
                    temp_frame=frame_hold{num2};
                    val_max(num2)=temp_frame.max_disp;
                end
                [obj.max_disp,obj.max_disp_frame]=max(val_max);%maximum displacement
                for num3=1:frames %minimum displacement
                    temp_frame=frame_hold{num3};
                    val_min(num3)=temp_frame.min_disp;
                end
                [obj.min_disp,obj.min_disp_frame]=min(val_min);%minimum displacement
                for num4=1:frames %maximum velocity
                    temp_frame=frame_hold{num4};
                    val_vel(num4)=temp_frame.max_vel;
                end
                %store values
                [obj.max_velocity,obj.max_velocity_frame]=max(val_vel);
                obj.num_frame=frames;
                obj.scale=scale;
                obj.timestep=timestep;
                obj.frame_holder=frame_hold;
                cd(directory) %return to directory
            end
        end %constructor
        %plotting functions
        function contour_frame=plot_contour(obj,select_string,ncont,mapstyle)
        %PLOT_CONTOUR Creates contour plot of displacement or velocity from vidobj. 
        %
        %Usage: contour_frame=plot_contour(obj,select_string,ncont,mapstyle)
        %
        %Inputs: select_string (string) specifies plotting of displacement
        %'disp' or velocity 'vel' data. ncont (scalar) specifies the number
        %of contours to use in the plot. mapstyle specifies the colormap to
        %plot with. Any of the Matlab color maps may be specified.
        %
        %Returns: a 1xn Matlab movie frame structure with fields cdata and
        %colormap.
        %
        %See also:
        %PLOT_QUIVER
        
            figure
            switch select_string
                case 'disp'%displacement magnitude contour
                    for n=1:obj.num_frame
                        contourf(obj.x_mat,obj.y_mat,obj.frame_holder{n}.disp_mat,ncont,'LineStyle','none');%plot displacement magnitude
                        axis ij
                        axis([0 obj.xdim 0 obj.ydim])
                        caxis([0,obj.max_disp])
                        colormap(mapstyle)
                        colorbar;
                        title(['Contraction displacement (\mum)' ' Frame: ' num2str(n)]);
                        contour_frame(n)=getframe(gcf);
                    end
                case 'vel' %velocity magnitude contour
                    for n=1:obj.num_frame
                        contourf(obj.x_mat,obj.y_mat,obj.frame_holder{n}.vel_mat,ncont,'LineStyle','none');%plot displacement magnitude
                        axis ij
                        axis([0 obj.xdim 0 obj.ydim])
                        caxis([0,obj.max_velocity])
                        colormap(mapstyle)
                        colorbar;
                        title(['Contraction velocity (\mum/s)' ' Frame: ' num2str(n)]);
                        xlabel('x');
                        ylabel('y');
                        contour_frame(n)=getframe(gcf);
                    end
            end
        end %plot_contour
        function quiver_frame=plot_quiver(obj,select_string,qscale,mapstyle,cscale)
        %PLOT_QUIVER Creates quiver plot of displacement or velocity from vidobj. 
        %
        %Usage: quiver_frame=plot_quiver(obj,select_string,qscale,mapstyle,cscale)
        %
        %Inputs: select_string (string) specifies plotting of displacement
        %'disp' or velocity 'vel' data. qscale (scalar) specifies the scale
        %of the quiver arrows in the plot. mapstyle specifies the colormap
        %to plot with. Any of the Matlab color maps may be specified.
        %cscale should either be string 'auto' for autoscaling the color
        %axis or vector [cmin,cmax] to define the min and max range of the
        %color axis.
        %
        %Returns: a 1xn Matlab movie frame structure with fields cdata and
        %colormap.
        %
        %See also:
        %PLOT_CONTOUR
        %QUIVER_MOD
        
            figure
            switch select_string
                case 'disp'%displacement quiver plot
                    for n2=1:obj.num_frame
                        if ischar(cscale) %maximum displacement as color scale
                            mincaxis=0;
                            maxcaxis=obj.max_disp;
                        else %use specified limits for color scale
                            mincaxis=cscale(1);
                            maxcaxis=cscale(2);
                        end
                        q=quiver_mod(obj.x_mat,obj.y_mat,obj.frame_holder{n2}.xdisp,obj.frame_holder{n2}.ydisp,qscale,mapstyle,mincaxis,maxcaxis);
                        set(gca,'Color','k')% for clear: set(gca,'color','none')
                        title(['Displacement [\mum]' ' Frame: ' num2str(n2)]);
                        caxis([mincaxis maxcaxis])
                        axis([0 obj.xdim 0 obj.ydim])
                        colorbar;
                        quiver_frame(n2)=getframe(gcf);
                    end
                case 'vel' %velocity quiver plot
                    for n2=1:obj.num_frame
                        if ischar(cscale)
                            mincaxis=0;
                            maxcaxis=obj.max_velocity;
                        else
                            mincaxis=cscale(1);
                            maxcaxis=cscale(2);
                        end
                        q=quiver_mod(obj.x_mat,obj.y_mat,obj.frame_holder{n2}.u_mat,obj.frame_holder{n2}.v_mat,qscale,mapstyle,mincaxis,maxcaxis);
                        set(gca,'Color','none')% for clear: set(gca,'color','none')
                        title(['Velocity (\mum/s)' ' Frame: ' num2str(n2)]);
                        caxis([mincaxis maxcaxis])
                        axis([0 obj.xdim 0 obj.ydim])
                        colorbar;
                        quiver_frame(n2)=getframe(gcf);
                    end
            end
        end %plot_quiver
        function threshold_frame=plot_threshold_contour(obj,ncont,mapstyle)
        %PLOT_THRESHOLD_CONTOUR Creates contour plot of displacement values above threshold. 
        %Prior to plotting the function mask_threshold must first be used
        %on the vidobj to create a binary mask.
        %
        %Usage: threshold_frame=plot_threshold_contour(obj,ncont,mapstyle)
        %
        %Input: ncont (scalar) specifies the number of contours to use in
        %the plot. mapstyle specifies the color map to plot with. Any of
        %the Matlab colormaps may be specified.
        %
        %Returns: threshold_frame, a 1xn Matlab movie frame structure with
        %fields cdata and colormap.
        
            figure %threshold plot
            for n3=1:obj.num_frame
                temp_framearea=obj.frame_holder{n3}.disp_mat;
                test=obj.frame_holder{n3}.binary_mask;
                filtered=test.*temp_framearea;
                h=contourf(obj.x_mat,obj.y_mat,filtered,ncont,'LineStyle','none');
                axis ij
                axis([0 obj.xdim 0 obj.ydim])
                caxis([0,obj.max_disp])
                colormap(mapstyle)
                colorbar;
                title(['Above threshold'  ' Frame: ' num2str(n3)]);
                threshold_frame(n3)=getframe(gcf);
            end
        end %plot_threshold_contour
        function figure_handle=heatmap_frequency(obj,ncont,mapstyle)
        %HEATMAP_FREQUENCY Creates heatmap of contraction frequency from vidobj.
        %The thresholded binary masks are summed togetherr and normalized
        %by the number of frames.
        %
        %Usage: figure_handle=heatmap_frequency(obj,ncont,mapstyle)
        %
        %Inputs: ncont (scalar) specifies the number of contours to use in
        %the plot. mapstyle specifies the color map to use. Any of the
        %Matlab color maps may be used.
        %
        %Returns: figure handle.
        %
        %See also:
        %HEATMAP_MAGNITUDE
            
            mask=zeros(size(obj.frame_holder{1,1}.binary_mask)); %initialize summed mask
            for n=1:length(obj.frame_holder)%add all masks together
                mask=mask+obj.frame_holder{1,n}.binary_mask;
            end
            mask=mask./n;%normalize mask
            figure_handle=figure;
            contourf(obj.x_mat,obj.y_mat,mask,ncont,'LineStyle','none')
            axis([0 obj.xdim 0 obj.ydim])
            ax=gca;
            hold on
            axis ij
            set(ax,'LineWidth',2)
            set(ax,'FontSize',30)
            colormap(mapstyle)
            colorbar
            title('Frequency heat map')
        end %heatmap_contracting
        function figure_handle=heatmap_magnitude(obj,ncont,mapstyle)
        %HEATMAP_MAGNITUDE Creates heatmap of contraction magnitude from vidobj. 
        %The displacement matrices are summed together and normalized by
        %the number of frames.
        %
        %Usage: figure_handle=heatmap_magnitude(obj,ncont,mapstyle)
        %
        %Inputs: ncont (scalar) specifies the number of contours to use in
        %the plot. mapstyle specifies the color map to use. Any of the
        %Matlab color maps may be used.
        %
        %Returns: figure handle.
        %
        %See also:
        %HEATMAP_FREQUENCY
        
            mask=zeros(size(obj.frame_holder{1,1}.disp_mat)); %initialize
            for n=1:length(obj.frame_holder)%Sum all frames
                mask=mask+obj.frame_holder{1,n}.disp_mat;
            end
            mask=mask./n;%normalize
            figure_handle=figure;
            contourf(obj.x_mat,obj.y_mat,mask,ncont,'LineStyle','none')
            axis([0 obj.xdim 0 obj.ydim])
            ax=gca;
            hold on
            axis ij
            set(ax,'LineWidth',2)
            set(ax,'FontSize',30)        
            colormap(mapstyle)
            colorbar
            title('Magnitude heat map')
        end %contraction_mag_heatmap
        function surf_frame=surfplot(obj,mapstyle)
        %SURFPLOT Visualizes displacement in 3D plot using z axis to encode magnitude
        %
        %Usage: surf_frame=surfplot(obj,mapstyle)
        %
        %Inputs: mapstyle specifies color map. Any of the Matlab color maps
        %may be used.
        %
        %Returns: a 1xn Matlab movie frame structure with fields cdata and
        %color map.
        
            AZ=-37.5;%default: -37.5
            EL=30; %default: 30
            figure
            for n2=1:obj.num_frame
                %create 3D plot
                h=surf(obj.x_mat,obj.y_mat,obj.frame_holder{n2}.disp_mat);
                %formatting here:
                title(['Displacement [\mum]' ' Frame: ' num2str(n2)]);
                axis([0 obj.xdim 0 obj.ydim 0 obj.max_disp])
                axis ij
                xlabel('x');
                ylabel('y');
                zlabel('z');
                caxis([0,obj.max_disp])
                colormap(mapstyle)
                colorbar;
                %adjust appearance and view
                material shiny
                view(AZ,EL)
                %view(AZ+n2*2,EL) %uncomment to rotate
                shading flat
                lighting flat
                lightangle(-45,30)
                h.AmbientStrength=0.8;
                h.DiffuseStrength=0.8;
                h.SpecularStrength=0.9;
                h.SpecularExponent=40;
                h.BackFaceLighting='unlit';
                surf_frame(n2)=getframe(gcf);
            end
        end %surfplot
        function fig_hand=plot_contraction_volume(obj)
        %PLOT_CONTRACTION_VOLUME Plots the contraction volume against time for the vidobj. 
        %
        %Usage: fig_hand=plot_contraction_volume(obj)
        %
        %Returns: figure handle.
        %
        %See also:
        %CONTRACTION_VOLUME
        
            frames=obj.num_frame;
            volumes=zeros(1,frames);%preallocate
            time=0:obj.timestep:(frames-1)*obj.timestep; %time vector
            for n=1:frames %get contraction volume for each frame
                volumes(1,n)=obj.frame_holder{1,n}.contraction_volume;
            end
            fig_hand=figure;
            plot(time,volumes)
            title('Contraction volume')
            ylabel('(\mum^3)')
            xlabel('time (s)')
        end %plot_contraction_volume
        function video=plot_binary_mask(obj)
        %PLOT_BINARY_MASK Animation formed by plotting each thresholded binary mask of vidobj.
        %Thresholded contraction is plotted as white on a black background.
        %
        %Usage: video=plot_binary_mask(obj)
        %
        %Returns: a 1xn Matlab movie frame structure with fields cdata and
        %color map.
        
            figure %threshold plot
            for n3=1:obj.num_frame %get each binary mask and plot
                binary_frame=obj.frame_holder{1,n3}.binary_mask;
                imshow(binary_frame,'InitialMagnification',600);%plot velocity magnitude
                axis xy
                video(n3)=getframe(gcf);
            end
        end %plot_binary_mask
        function video=xy_quiver(obj,qscale)
        %XY_QUIVER Animates quiver plot of displacement components.
        %Displacement in the x direction is plotted with blue
        %arrows and displacement in y direction with red arrows.
        %
        %Usage: video=xy_quiver(obj,qscale)
        %
        %Inputs: qscale (scalar) specifies the quiver arrow scale. 
        %
        %Returns a 1xn Matlab movie frame structure with fields cdata and
        %color map.
        %
        %See also:
        %PLOT_QUIVER
        
            zero_mat=zeros(size(obj.x_mat));%initialize zeros matrix
            figure
            for n=1:obj.num_frame %plot x and y displacement in quiver
                quiver(obj.x_mat,obj.y_mat,obj.frame_holder{n}.xdisp,zero_mat,qscale,'b');
                axis ij
                set(gca,'Color','k')
                hold on %overlay
                quiver(obj.x_mat,obj.y_mat,zero_mat,obj.frame_holder{n}.ydisp,qscale,'r');
                hold off
                title(['X disp. blue, Y disp. red [\mum]' ' Frame: ' num2str(n)]);
                axis([0 obj.xdim 0 obj.ydim])
                video(n)=getframe(gcf);
            end
        end %xy_quiver
        function fig_hand=plot_centroid_areathreshold(obj,area_threshold)
        %PLOT_CENTROID_AREATHRESHOLD Plots centroids for areas greater than area_threshold (scalar). 
        %This helps to check the spread of the larger contracting areas.
        %
        %Usage: fig_hand=plot_centroid_areathreshold(obj,area_threshold)
        %
        %Inputs: area_threshold (scalar): areas above this threshold are
        %included in the plot.
        %
        %Returns: figure handle.
        
            measure_prop='Centroid';
            [yind,~]=size(obj.x_mat);
            fig_hand=figure;
            hold on
            for fr=1:length(obj.frame_holder)%get morphology prop from each frame
                frame1=obj.frame_holder{1,fr}.morphology_props;
                for n=1:length(frame1)
                    area_measure=frame1(n).Area;%get area for each frame
                    if area_measure>area_threshold %check if area is greater than threshold
                        cent=frame1(n).(measure_prop);%get centroid
                        cent(2)=yind-cent(2);
                        cent=cent.*obj.vect_space.*obj.scale;
                        scatter(cent(1),cent(2),60,area_measure,'filled') %plot centroid
                    end
                end
            end
            hold off
            %plot formatting
            ax=gca;
            box(ax,'on')
            axis ij
            colorbar
            title('Centroids of thresholded area [\mum^2]')
        end %plot_centroid_areathreshold
        function fig_hand=plot_max_disp_each_frame(obj)
        %PLOT_MAX_DISP_EACH_FRAME Plots the maximum displacement in each frame.
        %Interpretation of this plot is limited since the trace is of the
        %maximum displacement and not of a single location in the frame.
        %Using the function plot_disp_location may be more beneficial.
        %
        %Usage: fig_hand=plot_max_disp_each_frame(obj)
        %
        %Returns: figure handle.
        %
        %See also:
        %PLOT_DISP_LOCATION
        
            val=zeros(1,length(obj.frame_holder));
            for n=1:length(obj.frame_holder) %plot max disp
                val(n)=obj.frame_holder{n}.max_disp;
            end
            t=(1:n)*obj.timestep;%time vector
            fig_hand=figure;%initialize figure with handle
            plot(t,val)
            xlabel('Time (s)')
            ylabel('Displacement (\mum)')
        end %plot_max_disp_each_frame
        function disp=plot_disp_location(obj,index)
        %PLOT_DISP_LOCATION Plots displacement for each frame at location specified by index.
        %
        %Usage: disp=plot_disp_location(obj,index)
        %
        %Input: index can either be string 'max' to plot at the global
        %maximum, 'min' to plot at the global minimum, or specify index as
        %vector [x,y] to plot at index location.
        %
        %Returns 1xn vector of displacement specified by displacement.
        
            if isnumeric(index) %if index is numeric, get row and column index
                row_ind=index(1);
                col_ind=index(2);
            else
                switch index
                    case 'max' %plot at overall maximum
                        row_ind=obj.frame_holder{1,obj.max_disp_frame}.disp_max_index(1);
                        col_ind=obj.frame_holder{1,obj.max_disp_frame}.disp_max_index(2);
                    case 'min' %plot at overall minimum
                        row_ind=obj.frame_holder{1,obj.min_disp_frame}.disp_min_index(1);
                        col_ind=obj.frame_holder{1,obj.min_disp_frame}.disp_min_index(2);
                end
            end
            disp=get_disp(obj,row_ind,col_ind); %call function to get displacement
            t=(1:length(disp))*obj.timestep;
            figure
            plot(t,disp)
            xlabel('Time (s)')
            ylabel('Displacement (\mum)')
        end %plot_disp_location
        
        %utility functions
        function area_vect=percent_area_disp_thresh(obj,threshold)
        %PERCENT_AREA_DISP_THRESH Obtains percent contracting area above displacement threshold
        %
        %Usage: area_vect=percent_area_disp_thresh(obj,threshold)
        %
        %Inputs: threshold (scalar), the value to use for the displacement
        %threshold.
        %
        %Returns: 1xn vector containing percent areas.
        
            for n=1:obj.num_frame%threshold each frame, calculate percent area
                temp_frame=obj.frame_holder{1,n};
                one_frame=temp_frame.disp_mat;%get displacement matrix
                binary_mask=one_frame>threshold;%threshold displacement
                area_vect(n)=100*(sum(binary_mask(:))/length(binary_mask(:)));
            end
        end %percent_area_disp_thresh
        function mask_threshold(obj,threshold)
        %MASK_THRESHOLD For each frame, obtains binary mask for contracting areas with displacements greater than threshold. 
        %Calls find_mask in frameobj which sets each binary mask frame.
        %
        %Usage: mask_threshold(obj,threshold)
        %
        %Inputs: threshold (scalar), the displacement threshold to use.
        %
        %See also:
        %FRAMEOBJ/FIND_MASK
        
            for n=1:obj.num_frame
                find_mask(obj.frame_holder{1,n},threshold)
            end
        end %mask_threshold
        function contraction_volume(obj)
        %CONTRACTION_VOLUME Loops through each frame and calculates the contraction volume. 
        %Calls frameobj method which sets contraction volume in the
        %frameobj properties.
        %
        %Usage: contraction_volume(obj)
        %
        %See also:
        %FRAMEOBJ/CALC_CONTRACTION_VOLUME
        
        
            for n=1:obj.num_frame
                calc_contraction_volume(obj.frame_holder{1,n});%frameobj method
            end
        end %contraction_volume
        function Area_props=measure_morphology(obj,type)
        %MEASURE_MORPHOLOGY Calculates morphology properties by passing each frame to calc_morphology, a frameobj method. 
        %The morphology properties are stored in the frameobj
        %"morphology_props" field. The function mask_threshold must be used
        %first to create binary masks which measure_morphology can then act
        %on.
        %
        %Usage: Area_props=measure_morphology(obj,type)
        %
        %Input: type (string) either 'box' or 'convex_hull'. Specifies if
        %regionprops should calculate BoundingBox or ConvexHull for later
        %use. If 'box' is specified then ConvexHull is not calculated. This
        %operation can be computationally expensive.
        %
        %Returns: cell array containing morphology data. 
        
            numFrames=obj.num_frame;
            B=cell(1,numFrames);%preallocate cells for holding bwconncomp output structure
            for n=1:numFrames %for each frame
                ind_frame=obj.frame_holder{1,n}; %get frame
                calc_morphology(ind_frame,obj.vect_space,type); %compute morphology
                B{n}=ind_frame.morphology_props;%set property in frameobj
            end
            Area_props=B;
        end%measure_morphology
        function calc_strain(obj,type) 
        %CALC_STRAIN Calculates strain tensor for all video frames. 
        %Passes each frame to frameobj method calc_tensor. The strain
        %values are stored in the frameobj field strain_tensor.
        %
        %Usage: calc_strain(obj,type)
        %
        %Inputs: type (string) either 'cauchy' or 'green_lagrangian' to use
        %the cauchy or green-lagrangian calculation methods respectively.
        %
        %See also:
        %FRAMEOBJ/CALC_TENSOR
        
            spacing=obj.x_mat(1,2)-obj.x_mat(1,1); %find spacing between points
            for n=1:obj.num_frame
                calc_tensor(obj.frame_holder{1,n},spacing,type)%a frameobj method
            end
        end%calc_strain
        function strain=transform_strain(obj,transform_type)
        %TRANSFORM_STRAIN Transforms the raw strain data using method specified by transform_type.
        %
        %Usage: strain=transform_strain(obj,transform_type)
        %
        %Inputs: transform_type (string), setting to 'principal' calculates
        %the principal strain from raw strain data. The function
        %calc_strain must be used first to create the raw strain data to
        %transform. This function to be expanded in future releases with
        %more transform options.
        %
        %Returns: 1xn cell array containing strain structure.
        %
        %See also:
        %CALC_STRAIN
        
            switch transform_type
                case 'principal'
                    strain=cell(1,obj.num_frame);
                    for n=1:obj.num_frame %calculate principal strain for each frame
                        strain{n}=calc_principal_strain(obj.frame_holder{1,n});
                    end
            end
        end %transform_strain
        function strain=strain_to_cell(obj)
        %STRAIN_TO_CELL Collects the strain tensor field from each frameobj and returns the strain structure in a 1xn cell array.
        %
        %Usage: strain=strain_to_cell(obj)
        %
        %Returns: strain, a 1xn cell array containing strain structure
        %
        %See also:
        %CALC_STRAIN
        
            strain=cell(1,obj.num_frame);%initialize cell array
            for n=1:obj.num_frame %get strain tensor from each frame
                strain{n}=obj.frame_holder{1,n}.strain_tensor;
            end
        end %strain_to_cell
        function [sorted_max,frame_index]=sort_morphology(obj,morph_prop)
        %SORT_MORPHOLOGY Sorts the morphology measurement specified by morph_prop from each frame from largest to smallest. 
        %
        %Usage: [sorted_max,frame_index]=sort_morphology(obj,morph_prop)
        %
        %Input: morph_prop (string) specifying the frameobj morphology
        %property to sort: 'Area', 'Perimeter', 'MajorAxisLength',
        %'MinorAxisLength','Eccentricity', 'Orientation',
        %'EquivDiameter','Centroid', 'percent_area'.
        %
        %Returns: sorted_max (1xn array) containing morphology values and
        %frame_index, a 1xn vector containing the index to the frame of the
        %morphology value.
        %
        %See also:
        %MEASURE_MORPHOLOGY
        
            frame_max=zeros(1,obj.num_frame); %preallocate
            switch morph_prop
                case 'percent_area'
                    for n=1:obj.num_frame %get percent_area from each frame
                        ind_frame=obj.frame_holder{1,n};
                        frame_max(n)=ind_frame.percent_area;
                    end
                otherwise
                    for n=1:obj.num_frame %iterate through frames to get property
                        ind_frame=obj.frame_holder{1,n};%get individual frame
                        frame_max(n)=get_max_morphology(ind_frame,morph_prop);
                    end
            end
            [sorted_max,frame_index]=sort(frame_max,'descend');
        end %sort_morphology
        function sorted=sort_all_areas(obj)
        %SORT_ALL_AREAS Sorts all areas in descending order.
        %Region_tree calls this method to sort areas into data tree.
        %
        %Usage: sorted=sort_all_areas(obj)
        %
        %Returns: sorted nx3 array with each row containing: area,
        %frame_index, morphology_index.
        %
        %See also:
        %REGION_TREE
        
            A=0;
            frame_index=0;
            morph_index=0;
            for n=1:obj.num_frame
                ind_frame=obj.frame_holder{1,n};
                [num_area,~]=size(ind_frame.morphology_props);
                while num_area>0
                    A(end+1)=ind_frame.morphology_props(num_area).Area;
                    frame_index(end+1)=n;
                    morph_index(end+1)=num_area;
                    num_area=num_area-1;
                end
            end
            A(1)=[];
            frame_index(1)=[];
            morph_index(1)=[];
            to_sort=[A;frame_index;morph_index]';
            [sorted,~]=sortrows(to_sort,-1); %sort rows by area column in descending order
            sorted=sorted';
        end %sort_all_areas
        function construct_tree(obj,type)
        %CONSTRUCT_TREE Passes vidobj to region_tree which constructs a data tree from the areas. 
        %Stores constructed tree in vidobj property tree. 
        %
        %Usage: construct_tree(obj,type)
        %
        %Input: type is a string either 'box' or 'convex_hull'. Specifies
        %method to use to check if child bounded by parent region. 'box'
        %may be significantly faster but may not handle complex shapes
        %well. 'convex_hull' is more accurate but may be significantly
        %slower. Future releases will add more options.
        %
        %See also:
        %REGION_TREE
        
            obj.tree=region_tree(obj,type);
        end
        function out=get_disp(obj,row_index,column_index)
        %GET_DISP Obtains displacement at index for every frame. 
        %
        %Usage: out=get_disp(obj,row_index,column_index)
        %
        %Inputs: row_index (scalar), the row index to query, column_index
        %(scalar), the column index to query.
        %
        %Returns: 1xn array containing displacement at location specified
        %by row_index and column_index.
        
            out=zeros(1,obj.num_frame);
            for n=1:obj.num_frame
                out(1,n)=obj.frame_holder{n}.disp_mat(row_index,column_index);
            end
        end
        function out=get_contraction_volume(obj) 
        %GET_CONTRACTION_VOLUME Gets contraction volume from all frames. 
        %
        %Usage: out=get_contraction_volume(obj)
        %
        %Returns 1xn array containing contraction volume for every frame.
        %
        %See also:
        %CONTRACTION_VOLUME
        
            out=zeros(1,obj.num_frame);%initialize
            for n=1:obj.num_frame%loop through each frame and get contraction volume
                out(n)=obj.frame_holder{n}.contraction_volume;
            end
        end  
        function out=get_max_velocity(obj)
        %GET_MAX_VELOCITY Gets maximum velocity from each frame.
        %
        %Usage: out=get_max_velocity(obj)
        %
        %Returns: 1xn array containing maximum velocity magnitude for every
        %frame.
        
            out=zeros(1,obj.num_frame);
            for n=1:obj.num_frame%loop through all frameobjs
                out(n)=obj.frame_holder{n}.max_vel;%obtain maximum velocity
            end
        end
    end %methods
end %classdef