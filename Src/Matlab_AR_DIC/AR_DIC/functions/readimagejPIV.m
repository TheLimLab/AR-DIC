function [x_mat,y_mat,xdisp_mat, ydisp_mat, disp_mat, u_mat, v_mat, vel_mat]...
    =readimagejPIV(varargin)
% READIMAGEJPIV Read data files from output of ImageJ PIV plugin.
%
% Usage: [x_mat,y_mat,xdisp_mat, ydisp_mat, disp_mat, u_mat, v_mat, vel_mat]...
% =readimagejPIV(filename, scale, timestep, xprev (optional), yprev (optional));
%
% Input: filename (string) the input file to read, scale (scalar) defines
% the pixel scale e.g. px/um, timestep (scalar) defines the time between
% frames, xprev (optional input, nxm array) contains the previous
% x-displacement matrix, yprev (optional input, nxm array) contains the
% previous y-displacement matrix.
%
% Returns: x_mat (nxm array) of x location values, y_mat (nxm array) of y
% location values, xdisp_mat (nxm array) of x displacements, ydisp_mat (nxm
% array) of y displacements, disp_mat nxm array of displacement magnitudes,
% u_mat: nxm array of velocity magnitude in x direction, v_mat: nxm array
% of velocity magnitude in y direction, vel_mat: nxm array of velocity
% magnitude
% 
% Format of txt files from ImageJ plugin PIV website manual:
% (https://sites.google.com/site/qingzongtseng/piv/tuto#post): x y ux1 uy1
% mag1 ang1 p1 ux2 uy2 mag2 ang2 p2 ux0 uy0 mag0 flag (x, y) is the
% position of the vector (center of the interrogation window). ux1, uy1 are
% the x and y component of the vector (displacement) obtained from the 1st
% correlation peak. mag1 is the magnitude (norm) of the vector. ang1, is
% the angle between the current vector and the vector interpolated from
% previous PIV iteration. p1 is the correlation value of the 1st peak. ux2,
% uy2, mag2, ang2, p2 are the values for the vector obtained from the 2nd
% correlation peak. ux0, uy0, mag0 are the vector value at (x, y)
% interpolated from previous PIV iteration. flag is a column used for mark
% whether this vector value is interpolated (marked as 999) or switched
% between 1st and 2nd peak (marked as 21), or invalid (-1).
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln

filename=varargin{1};
scale=varargin{2};
timestep=varargin{3};
M=dlmread(filename); %read in text file from ImageJ PIV
x_temp=M(:,1)./scale; %transfer data to individual vectors for coordinates
y_temp=M(:,2)./scale;
if nargin>3 %if displacement from previous frame is given use relative calculation
    xprev=varargin{4};
    yprev=varargin{5};
    xdisp_temp=M(:,3)./scale; %get x displacement and scale it by px/um ratio
    ydisp_temp=-M(:,4)./scale; %get y displacement and scale it. Flip sign for convention
    xdisp_temp=xdisp_temp+xprev;
    ydisp_temp=ydisp_temp+yprev;
else
xdisp_temp=M(:,3)./scale; %get x displacement and scale it by px/um ratio
ydisp_temp=-M(:,4)./scale; %get y displacement and scale it. Flip sign for convention
end
dispmagnitude_temp=sqrt((xdisp_temp.^2)+(ydisp_temp.^2)); %displacement magnitude
u_temp=xdisp_temp./timestep; %x velocity
v_temp=ydisp_temp./timestep; %y velocity
velmag_temp=sqrt((u_temp.^2)+(v_temp.^2));%calculate velocity magnitude

%reshape matrices for plotting and storage
firsty=y_temp(1);%find length of y dimension
n2=sum(y_temp==firsty);
n1=length(y_temp)/n2;%length of x dimension
x_mat=rot90(reshape(x_temp,n2,n1)); %rotate 90 degrees for plotting
y_mat=rot90(reshape(y_temp,n2,n1));
xdisp_mat=rot90(reshape(xdisp_temp,n2,n1));
ydisp_mat=rot90(reshape(ydisp_temp,n2,n1));
disp_mat=rot90(reshape(dispmagnitude_temp,n2,n1));
u_mat=rot90(reshape(u_temp,n2,n1));
v_mat=rot90(reshape(v_temp,n2,n1));
vel_mat=rot90(reshape(velmag_temp,n2,n1));
