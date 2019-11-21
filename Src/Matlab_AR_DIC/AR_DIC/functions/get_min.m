function [min_val,min_frame,row_ind,col_ind]=get_min(frames,min_type)
% GET_MIN Finds minimum strain value for type of strain specified by min_type.
%
% Usage: [min_val,min_frame,row_ind,col_ind]=get_min(frames,min_type)
%
% Input: frames is a 1xn cell array containing strain structure. min_type
% is string specifying which minimum to find.
%
% min_type options: 'principal_1' (minimum of primary principal strain),
% 'principal_2' (minimum principal strain), 'principal_min_shear' (minimum
% principal shear), 'raw_x' (minimum raw strain in x/horizontal direction,
% 'raw_y' (minimum raw strain in y/vertical direction), 'raw_xy' (minimum
% raw shear strain). The input frames can be obtained from the functions
% strain_to_cell or transform_strain which are both methods of the vidobj.
%
% Returns: minimum strain value min_val, for all frames including an index
% to the frame: min_frame, and the row and column indices row_ind, and
% col_ind.
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln
%
% See also:
% GET_MAX
% VIDOBJ/STRAIN_TO_CELL
% VIDOBJ/TRANSFORM_STRAIN

switch min_type
    case 'principal_1'%minimum of primary strain
        [r,c]=size(frames{1}.s1);
        for n=1:length(frames)
            [localmin(n),index(n)]=min(frames{n}.s1(:));
        end
        
    case 'principal_2' %minimum principal strain
        [r,c]=size(frames{1}.s2);
        for n=1:length(frames)
            [localmin(n),index(n)]=min(frames{n}.s2(:));
        end
        
    case 'principal_min_shear' %minimum transformed shear strain
        [r,c]=size(frames{1}.max_eng_shear);
        for n=1:length(frames)
            [localmin(n),index(n)]=min(frames{n}.max_eng_shear(:));
        end
        
    case 'raw_x' %minimum raw strain in x direction
        [r,c]=size(frames{1}.Exx);
        for n=1:length(frames)
            [localmin(n),index(n)]=min(frames{n}.Exx(:));
        end
        
    case 'raw_y' %minimum raw strain in y direction
        [r,c]=size(frames{1}.Eyy);
        for n=1:length(frames)
            [localmin(n),index(n)]=min(frames{n}.Eyy(:));
        end
    case 'raw_xy' %minimum raw shear strain
        [r,c]=size(frames{1}.Exy);
        for n=1:length(frames)
            [localmin(n),index(n)]=min(frames{n}.Exy(:));
        end
end
[min_val,min_frame]=min(localmin(:));
[row_ind,col_ind]=ind2sub([r,c],index(min_frame)); %get row, column indices for min_val
end