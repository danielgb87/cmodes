function dFC_utils_caps_create_nii(map_vector, mask_info, map_name)

% This function creates in the current folder, a .nii file for any given
% vector of activations, and mask. The vector length must be exactly equal
% to the amount of non-zero voxels in the mask.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   map_vector: stretched array of voxel activations.
%   mask_info:  structure with the mask binary map (as read by 
% OUTPUTS
%   data   : cell structure with a matrix (time x nodes(voxels)) for each subject with data normalized temporally.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 26-02-2020

%% 2. create a directory save maps
    
map = zeros(size(mask_info.mask));
ind = mask_info.mask_index;
for v = 1:length(map_vector)
    map(ind(1,v),ind(2,v),ind(3,v))  = map_vector(v);
end

V = spm_vol(mask_info.file);
V.fname = [pwd '/' map_name  '.nii'];
V.dt=[16,0];
spm_write_vol(V,map);
gzip([pwd '/' map_name  '.nii'])
delete ([pwd '/' map_name  '.nii'])

           

      

end