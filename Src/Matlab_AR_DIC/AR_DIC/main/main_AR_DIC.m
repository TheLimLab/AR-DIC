% Adaptive reference digital image correlation main script
% Script sections:
%   i. User input
%  ii. Adaptive Reference Digital Image Correlation (loop)
%       a. Digital image correlation
%       b. Adaptive reference evaluation
% iii. Reference score plotting
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% MAIN_AR_DIC_PROCESS
% MAIN_AR_DIC_PLOTTING
% MIJ_DIC
% DIC_disp

clear
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%User inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startframe=1; %set initial reference frame
endframe=160; %last video frame to analyze
%Match path to operating system '\' on Windows and '/' on Unix-based systems
folderpath_win='C:\Users\username\DIC_output';%set DIC output folder
folderpath='C:\\Users\\username\\DIC_output\\';%DIC output folder formatted for ImageJ %placing \\ after folder important
vidpath='C:\\Users\\username\\CM_contract.avi'; %video to analyze
method='fro'; %use frobenius norm method for adaptive reference comparison
norm_threshold=1.5;%frame becomes new reference if score is below norm threshold
stringtext='*PIV3*.txt'; %use only 3rd iteration
option_string='piv1=128 sw1=256 vs1=64 piv2=64 sw2=128 vs2=32 piv3=48 sw3=96 vs3=24 correlation=0.8';
%%%%%%%%%%%%%%%%%%%%%%%%%%End user inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%Adaptive reference digital image correlation%%%%%%%%%%%%%%%%%%%
%initialize, start MIJI
current_directory=pwd; %get current directory
currentframe=startframe+1;
refframe=startframe;
scorekeeper=[];
refstore=refframe;
Miji(false); %to start Fiji without gui: Miji(false);
%%
cd(folderpath_win) %change directory to get DIC output files
while currentframe<endframe+1 %Run AR-DIC while sufficient video frames exist
%DIC portion:
MIJ_DIC(folderpath,vidpath,refframe,currentframe,option_string);%perform DIC in ImageJ
%Adaptive Reference portion:
dispmag=DIC_disp(stringtext);%get DIC output from MIJI and format
refscore=calcrefscore(dispmag,method);%calculate reference score of current frame
scorekeeper(end+1)=refscore; %store reference score
%If reference score is below norm threshold then set current frame as reference
if refscore < norm_threshold
    refframe=currentframe; %set current frame as reference
    refstore(end+1)=refframe;%Store reference frame
end
currentframe=currentframe+1; %increment to next video frame
end
MIJ.exit; %Exit ImageJ
cd(current_directory); %return to original directory

%%
figure %plot reference score versus frame
plot(scorekeeper)%plot reference scores
hold on
scatter(refstore-1,zeros(1,length(refstore))); %mark reference frames
%EOF