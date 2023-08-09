function data = dFC_utils_postproc_liu_mask(data, mask_high,mask_low)

% This function masks frames by taking the highes and lowest percentiles of
% voxel intensities and zeroing the rest of voxels. This depending on a
% given high a low threshold percentile. This post-processing step improves
% signal-to-noise ration of the collected frames (ref [1]).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   data   : cell structure with a cell for every subject and data in
%   matrix form (time x nodes(voxels))
%   mask_high: percent of the highest voxels to choose
%   mask_low: percent of the lowest voxels to choose
% OUTPUTS
%   data   : cell structure with a matrix (time x nodes(voxels)) for each subject with data masked spatially.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%data = postproc_normalize(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[1] Liu, Xiao, and Jeff H. Duyn. "Time-varying functional network information 
% extracted from brief instances of spontaneous brain activity." Proceedings 
% of the National Academy of Sciences 110.11 (2013): 4392-4397.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 26-02-2020

Nvol = size(data,1);
p1 = mask_high; %% high threshold
p2 = mask_low;  %% low threshold
Nvox = size(data,2);

for i = 1:Nvol
    datavec = data(i,:);
    t1 = p1 * max(datavec); 
    t2 = p2 * min(datavec); 
    datavec(datavec<t1 & datavec>0) = 0;
    datavec(datavec>t2 & datavec<0) = 0;

    data(i,:) = datavec;
end

return;
