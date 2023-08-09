function dFC_prepare_data_parcel

% This function takes the input data, motion parameters, masks, and
% parcellation and extracts, preprocesses, and saves the data for dFC
% analysis.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 25-02-2020
% _________________________________________________________________________

clc,clear
%% EDIT HERE
% 1.a Define the data, parcellation, mask, and motion traces path, and
% their suffixes. Also remember that ROIs in the parcellation must by
% independent .nii files-
inputs.path_main = pwd;
inputs.path_data = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/13_smoothing';
inputs.path_mask = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/templates_masks/chd8_functional_template_mask.nii.gz';
inputs.path_parcel = '/home/safaai/DanielG/rois_atlas/rois_ludo/';
inputs.path_motion_data = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/07_motion_correction';   % leave equal to zero (inputs.path_motion_data = 0) if no motion analysis nor scrubbing is to be done

inputs.suffix_data= '*GqNegKNeg*_*pre_cno_bs*';   %to identify desired subjects
inputs.suffix_motion = '*GqNegKNeg*_*pre_cno_bs*_*.txt';    %to identify the desired motion traces.Leave equal to zero (inputs.suffix_motion = 0) if no motion analysis nor scrubbing is to be done
inputs.suffix_rois = '*.nii.gz';

% 1.b Define which postprocessing should be done on the data.
inputs.pproc_zscore = 1;    % Perform data normalizatuion (subtract mean and divide by STD)
inputs.pproc_scrub = 0;     % Perform data scrubbing based on framewise displacement flagged volumes.

% 1.c Define some basic data characteristics
inputs.analysis_TR = 1;   %temporal resolution (in seconds)
inputs.analysis_ID = 'gqneg_pre_cno_bs';    % a folder with this name will be created with the results

if isenum(inputs.path_motion_data) == 0
    % 1.d define the motion correction parameters.
        inputs.motion_radius           = 5;      % mouse brain radius (5mm by default)
        inputs.motion_divide           = 1;      % divide by 10 translation motion traces.
        inputs.motion_suite            = 'fsl'; % suite used to obtain motion parameters ('fsl', 'afni', or 'spm')
        inputs.motion_scrub_thr        = 0.1;   % Framewise displacement censoring criteria (in mm).
end
if inputs.pproc_scrub == 1
    % 1.e describe the scrubbing method
        inputs.motion_scrub_criteria    = 'FD';   % can be 'FD' or 'DVARS' depending on the criteria
        inputs.motion_scrub_plus    = [0];   % define if only scrub the flagged volume ([0]), or the previous and next volumes ([-1,0,1])    
end

%% 2. Extract data in matrix form within the selected mask, thextract the timecourses for each parcel.
% 2.1 extract data in matrix 
cd(inputs.path_main)
[data, mask_info, sub_list] = dFC_utils_prepare_dataset_matrix(inputs.path_data, inputs.suffix_data, inputs.path_mask);
inputs.Nsubs = length(sub_list);

% 2.2 parcellate data and ensure that no subject has an ROI with null data
cd(inputs.path_parcel)
inputs.parcel_roi_list = dir(fullfile(inputs.path_parcel,inputs.suffix_rois));
inputs.Nrois = length(inputs.parcel_roi_list);


for s = 1:length(roi_list)
    inputs.parcel_roi_id{s} = inputs.parcel_roi_list(s).name(1:end-7);
    roi_full_img = spm_read_vols(spm_vol(inputs.parcel_roi_list(s).name));

    for sub = 1:inputs.Nsubs
        vec = spm_vec(roi_full_img);
        vec = vec(spm_vec(mask_info.mask)>=1);
        roi_ind = find(vec>=1);
        data_parcel{sub,1}(:,s) = mean(data{sub}(:,roi_ind),2);
              
    end
end

%% 3. extract motion information
if isenum(inputs.path_motion_data) == 0
[motion_info] = dFC_utils_check_motion(inputs.path_motion_data, inputs.suffix_motion, inputs.motion_radius, inputs.motion_divide, inputs.motion_suite, inputs.motion_scrub_thr);
end

%% 4. scrub data.
if inputs.pproc_scrub == 1
data = postproc_censor_data(data, motion_info, inputs.motion_scrub_plus);
data_parcel = postproc_censor_data(data_parcel, motion_info, inputs.motion_scrub_plus);
end

%% 5. normalize data
if inputs.pproc_zscore==1
data = dFC_pproc_normalize_data(data);
data_parcel = dFC_pproc_normalize_data(data_parcel);
end

%% 6.save the inputs, motion aprameters, mask, and data in both formats (voxelwise and parcel)
cd(inputs.path_main)
save(['motion_info_' inputs.analysis_ID],'motion_info','-v7.3')
save(['inputs_' inputs.analysis_ID],'inputs')
save(['data_' inputs.analysis_ID],'data','-v7.3')
save(['data_parcel' inputs.analysis_ID],'data_parcel','-v7.3')





end %function