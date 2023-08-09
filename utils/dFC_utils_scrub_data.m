function data_censored = dFC_utils_scrub_data(data, motion_info, motion_scrub_plus) 

% This function scrubs data frames from the dataset using the motion_info
% structure from the motion check (see motion_check). It removes for each
% subject, frames that were flagged as being high motion depending on the
% criteria previously selected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   data   : cell structure with a cell for every subject and data in
%   matrix form (time x nodes(voxels))
%   motion_info: structure with the information for each subject regarding
%   motion. Specifically, the varialbe motion_info.FD_flag{subs} is used. 

% OUTPUTS
%   data_censored   : cell structure with a matrix (time x nodes(voxels)) for each subject prunned data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%data = postproc_normalize(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 25-02-2020

tmin = min(motion_scrub_plus);
tmax = max(motion_scrub_plus);
Nsubs = length(motion_info.FD_flag);
data_censored = cell(length(data),1);
flag = motion_info.FD_flag;

for sub = 1:Nsubs
    nobs = size(data{sub},1);
    
    data_censored{sub} = data{sub}(flag{sub}(1+abs(tmin):nobs-abs(tmax)) == 0,:);

end

end %function