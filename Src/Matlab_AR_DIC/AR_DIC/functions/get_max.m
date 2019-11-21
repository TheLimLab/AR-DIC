function [max_val,max_frame,row_ind,col_ind]=get_max(frames,max_type)
% GET_MAX Finds maximum strain value for type of strain specified by max_type.
%
% Usage: [max_val,max_frame,row_ind,col_ind]=get_max(frames,max_type)
%
% Input: frames is 1xn cell array containing strain structure. max_type is
% string specifying which maximum to find.
%
% max_type options: 'principal_1' (maximum principal strain), 'principal_2'
% (maximum of secondary principal strain), 'principal_max_shear' (maximum
% principal shear),'raw_x' (maximum strain in x/horizontal direction),
% 'raw_y' (maximum strain in y/vertical direction), 'raw_xy' (maximum raw
% shear strain). The input frames can be obtained from the functions
% strain_to_cell or transform_strain which are both methods of the vidobj.
%
% Returns: maximum strain value max_val, for all frames including an index
% to the frame: max_frame, and the row and column indices row_ind, and
% col_ind.
%
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% GET_MIN
% VIDOBJ/STRAIN_TO_CELL
% VIDOBJ/TRANSFORM_STRAIN

switch max_type
    case 'principal_1'%maximum principal strain
        [r,c]=size(frames{1}.s1);
        for n=1:length(frames)
            [localmax(n),index(n)]=max(frames{n}.s1(:));
        end
    case 'principal_2'%maximum of secondary strain
        [r,c]=size(frames{1}.s2);
        for n=1:length(frames)
            [localmax(n),index(n)]=max(frames{n}.s2(:));
        end
    case 'principal_max_shear'%maximum transformed shear strain
        [r,c]=size(frames{1}.max_eng_shear);
        for n=1:length(frames)
            [localmax(n),index(n)]=max(frames{n}.max_eng_shear(:));
        end
    case 'raw_x'%maximum raw strain in x direction
        [r,c]=size(frames{1}.Exx);
        for n=1:length(frames)
            [localmax(n),index(n)]=max(frames{n}.Exx(:));
        end
    case 'raw_y' %maximum raw strain in y direction
        [r,c]=size(frames{1}.Eyy);
        for n=1:length(frames)
            [localmax(n),index(n)]=max(frames{n}.Eyy(:));
        end
    case 'raw_xy' %maximum raw shear strain
        [r,c]=size(frames{1}.Exy);
        for n=1:length(frames)
            [localmax(n),index(n)]=max(frames{n}.Exy(:));
        end
end
[max_val,max_frame]=max(localmax(:));
[row_ind,col_ind]=ind2sub([r,c],index(max_frame)); %get row, column indices for max_val
end