function [data, mask_info, sub_list] = dFC_utils_prepare_dataset_matrix(data_path, pp_suffix, mask_path)


% This function organizes data from a dataset into a matrix form using a
% predefined mask. Be sure that the mask is in the same resolution space as
% the preprocessed data. Data in each matrix is determined by the amount of
% non-zero voxels in the input mask.

% Also make sure that data is organized as follows:
% data_path/subject_ID/preprocessing_step/subject_ID_suffix.nii
% if data is compressed, the program decompresses it, uses the uncompressed
% data then eliminates it automatically.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   data_path   : Path where subject data is stored. ex: 'D:\Gutierrez\mice_dFC\preprocessing_scan10_partA_500\subjects'
%   pp_suffix   : Suffix at the end of the data_file after the subjectID in case its not the same as the preproc_step. ex: 'smoothed'.
%   mask_path   : a binary .nii file path to an uncompressed mask. ex: 'D:\Gutierrez\mice_dFC\templates_for_Xan\template_hires_mask.nii'
% OUTPUTS
%   data        : cell structure with a matrix (time x nodes(voxels)) for each subject.
%   mask_info   : structure with relevant information about the mask including a map-version in single resolution and indexes of non-zero values.
%   sub_list    : cell with the ID of each subject.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Daniel Gutierrez-Barragan 2020, V1 25-02-2020

%% 1. load template mask information.

mask_file = mask_path;
mask = single(spm_read_vols(spm_vol(mask_file)));
% Check indexes with non-zero values
[mask_ind(1,:), mask_ind(2,:), mask_ind(3,:)] = ind2sub(size(mask),find(mask ==1));

mask_info.file = mask_file;
mask_info.mask = mask;
mask_info.mask_index = mask_ind;

%% 2. Load data and convert it to matrix form.
% Go to where subjects are and make a list of their IDs.
cd(data_path),
list = dir(pp_suffix);
Nsubs = length(list);

for sub = 1:Nsubs
    sub_list{sub,1} = list(sub).name;
end

% Load each subject and convert dta to matrix form. If data is zipped,
% unzip it, use the unzipped data and then eliminate it.
data = cell(Nsubs,1);

for sub=1:Nsubs
    
    % Go to the subjects preprocessing dir   
        temp = spm_vol(sub_list{sub,1});
        parfor t = 1:length(temp)
            data1(:,t) = spm_vec(spm_read_vols(temp(t)));
        end
        data{sub} = data1(spm_vec(mask)==1,:);
        data{sub,1} = data{sub,1}';

        clear data1 temp
            
end

% Deliver output data in single format (reduces 1/2 of memory usage).
for sub = 1:length(data)    
    data{sub,1} = single(data{sub});
end

end %function

