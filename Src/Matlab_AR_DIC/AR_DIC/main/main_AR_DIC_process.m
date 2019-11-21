% Process output of AR-DIC method: 'main_AR_DIC.m'
% 
% Organizes data structures, obtains mechanics measurements, and builds a
% data tree from the results. Data can be explored manually by viewing the
% vidobjs in the variable viewer. The script main_AR_DIC_plotting handles
% plotting of the results.
%
% Script sections:
%   i. User input
%  ii. Data processing
%       a. Initial processing in vidobjs
%       b. Secondary/specialized processing
% iii. Data exploration
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% MAIN_AR_DIC
% MAIN_AR_DIC_PLOTTING

clear
close all %close open figures
[~]=get(0,'Factory');%Commands to set all figure backgrounds to white
set(0,'defaultfigurecolor',[1 1 1]);
set(0,'defaultlinelinewidth',2);
clc

%%%%%%%%%%%%%%%%%%%%%% User Input: Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%
folder_stretch100='C:\Users\username\Documents\MATLAB\CM_contract'; %stretch100
folder_control100='C:\Users\username\Documents\MATLAB\CM_non_contract';
%-------------------------------------------------------------------------%
stringtext='*PIV3*.txt';% get only output files from 3rd iteration of plugin
scale100=1.04; %set pixel micrometer scale px/um (10x objective)
timestep=0.2; %set time between frames (seconds)
threshold=0.14; %define threshold for threshold plot 0.14 from Fukuda et. al
%plotting options:
ncont=50; %number of contour levels in plots
qscale=1; %quiver plot arrow scale
mapstyle='parula'; %map color scheme or try gray
strain_type='green_lagrangian'; % strain calculation to use. Options: 'green_lagrangian' or 'cauchy'
tree_opt='box';%'box' for faster processing, 'convex_hull' for higher accuracy
%%%%%%%%%%%%%%%%%%%%%%%%%%End of user input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Data process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Primary data processing occurs here:
%initial data processing occurs in vidobj, (displacement, velocity, etc.):
stretch2=vidobj(stringtext,folder_stretch100, scale100,timestep); %process contracting
control2=vidobj(stringtext,folder_control100, scale100,timestep); %process non-contracting
data={control2,stretch2};%combine into cell for iteration

%specialized data processing
for sets=1:length(data) %iterate for length of datasets
mask_threshold(data{sets},threshold) %set masks based on threshold
contraction_volume(data{sets})%calculate contraction volume
[~]=measure_morphology(data{sets},tree_opt); %calculate morphology parameters
construct_tree(data{sets},tree_opt); %build data tree of regions
[~]=data{sets}.tree.leaf_index; %get indices of all leaves (find smallest independent regions)
calc_strain(data{sets},strain_type);%calculate strains
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Explore data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User may wish to explore data here:
leaves_s100=get_node(stretch2.tree,stretch2.tree.leaf_index); %get nodes from leaf indices
leaf_path=trace_path(leaves_s100{1},'time_sort'); %trace path from specified leaf to root
[sorted_percent_area,percent_area_frame_index]=sort_morphology(stretch2,'percent_area');
[sorted_max,frame_index]=sort_morphology(stretch2,'Area');
%calculate principal strains:
strains_principal_s2=transform_strain(stretch2,'principal');
strains_principal_c2=transform_strain(control2,'principal');
strains_s2=strain_to_cell(stretch2);
strains_c2=strain_to_cell(control2);
%Find maximum, minimum strain
[max_val_s2,max_frame_s2,xind_s2,yind_s2]=get_max(strains_principal_s2,'principal_1');
[max_val_c2,max_frame_c2,xind_c2,yind_c2]=get_max(strains_principal_c2,'principal_1');
secondary_max_strain_s2=strains_principal_s2{max_frame_s2}.s2(xind_s2,yind_s2);
[min_val_s2,min_frame_s2,xind_min_s2,yind_min_s2]=get_min(strains_principal_s2,'principal_2');
[min_val_c2,min_frame_c2,xind_min_c2,yind_min_c2]=get_min(strains_principal_c2,'principal_2');
secondary_min_strain_s2=strains_principal_s2{min_frame_s2}.s1(xind_min_s2,yind_min_s2);

%uncomment to check raw strain max:
% [max_val_raw_s2,max_frame_raw_s2,xind_raw_s2,yind_raw_s2]=get_max(strains_s2,'raw_y');
% [min_val_raw_s2,min_frame_raw_s2,xind_min_raw_s2,yind_min_raw_s2]=get_min(strains_s2,'raw_y');
% [max_val_c2,max_frame_c2,xind_c2,yind_c2]=get_max(strains_c2,'raw_x');

%Run plotting file:
% run main_AR_DIC_plotting.m
%EOF