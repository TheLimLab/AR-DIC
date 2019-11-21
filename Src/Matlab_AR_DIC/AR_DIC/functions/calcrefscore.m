function refscore=calcrefscore(disp_mat,refmethod)
% CALCREFSCORE Calculates reference frame score from displacement matrix.
%
% Usage: refscore=calcrefscore(disp_mat,refmethod)
%
% Inputs: disp_mat, nxm displacement array. refmethod (string) specifies
% the scoring method. 'fro' uses the frobenius norm and 'ent' uses the
% entropy.
%
% Returns: score value refscore (scalar).
% 
% Adaptive Reference Digital Image Correlation v 1.0 2018
% Biomaterials and Mechanotransduction Lab University of Nebraska-Lincoln

switch refmethod
    case 'fro'
       refscore=norm(disp_mat,'fro');%frobenius norm
    case 'ent'
       refscore=entropy(disp_mat); %entropy
end