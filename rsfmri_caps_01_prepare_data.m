function caps_01_dFC_prepare_data

% This function takes the input data, motion parameters, masks, and
% parcellation and extracts, preprocesses, and saves the data for dFC
% analysis. It performs it on independent datasets. IF SCRUBBING IS
% PERFORMED, THE FUNCTION YIELDS THE TWO VERSIONS OF THE DATASET WITH THE
% SUFIX "_scrub_"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V2 27-01-2021
% _________________________________________________________________________

clc,clear
%% EDIT HERE
% path to CAPs_scripts and external toolboxes
addpath(genpath('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127'))

% 1.a Define the data, parcellation, mask, and motion traces path, and
% their suffixes. Also remember that ROIs in the parcellation must by
% independent .nii files-
inputs.path_main = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221';
addpath(genpath(inputs.path_main))

ndatasets = 2;
inputs.path_mask = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/chd8_functional_template_mask_woCereb_Vent_dgb.nii.gz';
inputs.path_parcel = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/atlas_by_networks/rois';

%1.b for each dataset, edit the following
inputs.path_motion_data{1} = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/example_data_awake/motion_traces';   % leave equal to zero (inputs.path_motion_data = 0) if no motion analysis nor scrubbing is to be done
inputs.path_data{1} = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/example_data_awake/smoothing';
inputs.suffix_data{1}= '*smoothed.nii*';   %to identify desired subjects
inputs.ds_name{1}='AWK';        %name the dataset

inputs.path_motion_data{2} = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/example_data_xan/motion_traces';   % leave equal to zero (inputs.path_motion_data = 0) if no motion analysis nor scrubbing is to be done
inputs.path_data{2} = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/example_data_xan/smoothing';
inputs.suffix_data{2}= '*smoothed.nii*';   %to identify desired subjects
inputs.ds_name{2}='ANE';


inputs.suffix_motion = '*mcf.txt';    %to identify the desired motion traces.Leave equal to zero (inputs.suffix_motion = 0) if no motion analysis nor scrubbing is to be done
inputs.suffix_rois = '*.nii*';

% 1.c Define which postprocessing should be done on the data. It is
% suggested to scrub the data (inputs.pproc_scrub=1) in order to prepare it
% if needed for posterior analyses. If you scrub, the unscrubbed data
% output will still be delivered
inputs.pproc_zscore = 1;    % Perform data normalizatuion (subtract mean and divide by STD)
inputs.pproc_scrub = 1;     % Perform data scrubbing based on framewise displacement flagged volumes.

% 1.d Define some basic data characteristics for each dataset and the
% motion correcton strategy used
inputs.analysis_TR = [1, 1.2];   %temporal resolution (in seconds)
inputs.analysis_ID = 'awk_ane';    % a folder with this name will be created with the results


if isenum(inputs.path_motion_data) == 0
    % 1.e define the motion correction parameters.
    inputs.motion_radius           = 5;      % mouse brain radius (5mm by default)
    inputs.motion_divide           = 1;      % divide by 10 translation motion traces.
    inputs.motion_suite            = 'fsl'; % suite used to obtain motion parameters ('fsl', 'afni', or 'spm')
    inputs.motion_scrub_thr        = 0.0750;   % Framewise displacement censoring criteria (in mm).
end
if inputs.pproc_scrub == 1
    % 1.f describe the scrubbing method
    inputs.motion_scrub_criteria    = 'FD';   % can be 'FD' or 'DVARS' depending on the criteria
    inputs.motion_scrub_plus    = [0];   % define if only scrub the flagged volume ([0]), or the previous and next volumes ([-1,0,1])
end

%% 2. Extract data in matrix form within the selected mask, thextract the timecourses for each parcel.
% 2.a extract data in matrix
cd(inputs.path_main)
for i = 1:ndatasets
    [data_vw{i}, mask_info, sub_list{i}] = dFC_utils_prepare_dataset_matrix(inputs.path_data{i}, inputs.suffix_data{i}, inputs.path_mask);
    inputs.Nsubs{i} = length(sub_list{i});
end

% 2.b parcellate data and ensure that no subject has an ROI with null data
cd(inputs.path_parcel)
inputs.parcel_roi_list = dir(fullfile(inputs.path_parcel,inputs.suffix_rois));
inputs.Nrois = length(inputs.parcel_roi_list);

inputs.flag_roisNAN = zeros(inputs.Nrois,ndatasets);
for nd = 1:ndatasets
    inputs.flag_roisNAN_sub{nd} = zeros(inputs.Nsubs{nd},length(inputs.parcel_roi_list));
    for s = 1:inputs.Nrois
        inputs.parcel_roi_id{s} = inputs.parcel_roi_list(s).name(1:end-7);
        roi_full_img = spm_read_vols(spm_vol(inputs.parcel_roi_list(s).name));
        inputs.roi_nvox(s) = sum(spm_vec(roi_full_img));
        vec = spm_vec(roi_full_img);
        vec = vec(spm_vec(mask_info.mask)>=1);
        roi_ind = find(vec>=1);
        inputs.roi_nvox_inmask(s) = length(roi_ind);
        for sub = 1:inputs.Nsubs{nd}
            
            data_parcel{nd}{sub,1}(:,s) = mean(data_vw{nd}{sub}(:,roi_ind),2);
            data_parcel{nd}{sub,1}(isnan(data_parcel{nd}{sub,1}))=0;
            Nobs = length(data_parcel{nd}{sub,1}(:,s));
            if sum(data_parcel{nd}{sub,1}(:,s))==0
                inputs.flag_roisNAN_sub{nd}(sub,s)=1;
            end
            data_parcel{nd}{sub,1}(:,s)=squeeze(data_parcel{nd}{sub,1}(:,s));
        end
        if sum(inputs.flag_roisNAN_sub{nd}(:,s))>0
            inputs.flag_roisNAN(s,nd)=1;
        end
    end
end


%% 3. extract motion information
if isenum(inputs.path_motion_data) == 0
    for i = 1:ndatasets
        [motion_info{i}] = dFC_utils_check_motion(inputs.path_motion_data{i}, inputs.suffix_motion, inputs.motion_radius, inputs.motion_divide, inputs.motion_suite, inputs.motion_scrub_thr);
    end
end

%% 4. scrub data.
if inputs.pproc_scrub == 1
    for nd = 1:ndatasets
        data_vw_scrub{nd} = dFC_utils_scrub_data(data_vw{nd}, motion_info{nd}, inputs.motion_scrub_plus);
        data_parcel_scrub{nd} = dFC_utils_scrub_data(data_parcel{nd}, motion_info{nd}, inputs.motion_scrub_plus);
    end
end

%% 5. normalize data
if inputs.pproc_zscore==1
    for nd = 1:ndatasets
        data_vw{nd} = dFC_pproc_normalize_data(data_vw{nd});
        data_parcel{nd} = dFC_pproc_normalize_data(data_parcel{nd});
    end
    if inputs.pproc_scrub==1
        for nd = 1:ndatasets
            data_vw_scrub{nd} = dFC_pproc_normalize_data(data_vw_scrub{nd});
            data_parcel_scrub{nd} = dFC_pproc_normalize_data(data_parcel_scrub{nd});
        end
    end
end

%% 6.save the inputs, motion aprameters, mask, and data in both formats (voxelwise and parcel)
cd(inputs.path_main)
save(['motion_info_0750_' inputs.analysis_ID],'motion_info','-v7.3')
save(['inputs_' inputs.analysis_ID],'inputs')
save(['data_vw_' inputs.analysis_ID],'data_vw','-v7.3')
save(['data_parcel_' inputs.analysis_ID],'data_parcel','-v7.3')
save(['mask_info_' inputs.analysis_ID],'mask_info','-v7.3')
save(['subject_list_' inputs.analysis_ID],'sub_list','-v7.3')
if inputs.pproc_scrub ==1
    save(['data_vw_scrub_' inputs.analysis_ID],'data_vw_scrub','-v7.3')
    save(['data_parcel_scrub_' inputs.analysis_ID],'data_parcel_scrub','-v7.3')
end




end %function
