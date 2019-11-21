function dispmag=DIC_disp(stringtext)
% DIC_DISP Calculates displacement magnitude from highest numerical DIC filename in directory.
%
% Usage: dispmag=DIC_disp(stringtext)
%
% Inputs: stringtext (string) specifies the text to search for to input DIC
% files.
%
% Returns: dispmag (array) containing displacement values from most recent
% DIC iteration.
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% MIJ_DIC
% VIDOBJ
% READIMAGEJPIV

PIVtext=dir(stringtext);%get text file list
PIVtext_names={PIVtext.name};%get file names
[~,idx]=sort_nat(PIVtext_names);%get sorting index
PIVtext=PIVtext(idx);%sort files in alphanumerical order
frames=length(PIVtext);%get total number of frames
filename=PIVtext(frames).name;%get most recent frame
M=dlmread(filename); %read in text file from ImageJ PIV
dispmag=M(:,5); %transfer data to individual vectors for coordinates