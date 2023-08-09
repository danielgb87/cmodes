function caps_utils_create_nii_parcel(map_vector, roi_info, mask_info, map_name)

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

% USAGE EXAMPLE:
% caps_create_nii(cap1, mask_info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 14-11-17

%% 2. create a directory save maps
    
map = zeros(size(mask_info.mask));
ind = mask_info.mask_index;

for r = 1:length(roi_info.vox_ind)
    for v = 1:size(roi_info.vox_ind{r},2)
        map(roi_info.vox_ind{r}(1,v),roi_info.vox_ind{r}(2,v),roi_info.vox_ind{r}(3,v))  = map_vector(r);
    end
end

V = spm_vol(mask_info.file);
V.fname = [pwd '/' map_name  '.nii'];
V.dt=[16,0];
spm_write_vol(V,map);

           

      

end