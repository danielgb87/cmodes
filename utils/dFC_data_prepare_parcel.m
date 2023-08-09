function dFC_data_prepare_parcel

% This function takes the input data, motion parameters, masks, and
% parcellation and extracts, preprocesses, and saves the data for dFC
% analysis.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 25-02-2020
% _________________________________________________________________________


%% EDIT HERE
% 1.a Define the data, parcellation, mask, and motion traces path, and
% their prefixes.
inputs.main = pwd;
data_path = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/13_smoothing';
mask_path = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/templates_masks/chd8_functional_template_mask.nii.gz';
motion_preproc = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/07_motion_correction';

data_suffix = '*GqNegKNeg*_*pre_cno_bs*'; %to identify desired subjects
motion_suffix = '*GqNegKNeg*_*pre_cno_bs*_*.txt'; %to identify the desired motion traces.



cd(inputs.main)
inputs.TR=1;
%% 1. Prepare data for analysis and check_motion (Takes around 0.5 minutes per subject with 500 timepoints and around 80.000 in-brain voxels).

% define which dataset are you using and its path; the preprocessing step
% and suffix; and the mask and its path.
data_path = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/';
inputs.data_path = data_path;
mask_path = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/templates_masks/chd8_functional_template_mask.nii.gz';
inputs.mask_path = mask_path;
pp_step   = '13_smoothing';
inputs.pp_step = pp_step;
pp_suffix = '*GqNegKNeg*_*pre_cno_bs*';
inputs.pp_suffix = pp_suffix;

%extract data into matrix form.
cd(inputs.main)
[data, mask_info, subject_list] = dataset_matrix_convert_in_folder(data_path, pp_step, pp_suffix, mask_path);


cd(inputs.main)
dataset_name = 'data_gqNeg_pre_cno';
inputs.dataset_name = dataset_name;
mask_info_name = 'mask_low_res_wo_Cereb_Vent';
inputs.mask_info_name = mask_info_name;
motion_info_name = 'motion_info_100_gqNeg_pre_cno';
inputs.motion_info_name = motion_info_name;
subject_list_name = 'subject_list_gqNeg_pre_cno';
inputs.subject_list_name = subject_list_name;


save(dataset_name,'data','-v7.3')
save(mask_info_name,'mask_info','-v7.3')
save(subject_list_name,'subject_list','-v7.3')



cd(inputs.main)
disp('done preparing data')

% Motion check data (takes around 3 second). Build an inputs structure with
% the motion checking parameters.
inputs.motion.radius           = 5;      % mouse brain radius (5mm by default)
inputs.motion.divide           = 1;      % divide by 10 translation motion traces.
inputs.motion.suite            = 'fsl'; % suite used to obtain motion parameters ('fsl', 'afni', or 'spm')
inputs.motion.censor_method    = 'FD';   % can be 'FD' or 'DVARS' depending on the criteria
inputs.motion.scrub_thr        = 0.05;   % Framewise displacement censoring criteria (in mm).
inputs.motion.DVARSscrub_thr   = 0.5;    % DVARS censoring criteria (in % signal change).
% define preprocessing step where the motion parameters are.
motion_preproc = '/media/DATA/DanielG/dbh_KORD_tonini/Preprocessing/07_motion_correction/';
motion_suffix  = '*GqNegKNeg*_*pre_cno_bs*_*.txt';

motion_info = motion_check_in_folder(inputs.motion, data_path,motion_preproc,motion_suffix);
inputs.pproc.motion_check = 1;
motion_info.inputs = inputs.motion;

disp('done checking motion')
cd(inputs.main)
save(motion_info_name,'motion_info','-v7.3')
save('inputs_gqNeg_pre_cno','inputs')

% another option is to load data already in cell-matrix format and a
% mask_info file. Remember that the mask_info structure MUST contain at
% least: 1-The mask_file path; 2- the 3D binary tensor (as given by
% spm_readvol); and 3- the three arrays of the indexes of non-zero
% elements.
%
load(inputs.dataset_name);
load(inputs.mask_info_name);
load(inputs.motion_info_name);
load(inputs.subject_list_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Post-processing of data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each step done, try to also save an inputs structure to save at the
% end of the analysis to have record of all processes done.

Nsubs = length(data);


%2.0 If only a part of the scan will be used, partition the data.
% for sub = 1:length(data)
%     data1{sub,1} = data{sub,1}(501:1000,:);
% end
% clear data
% data = data1; clear data1
% save('data_xan_xan_hires_pt2','data','-v7.3')

%2.0.1 if you want to notch filter the data
% Nsubs = length(data);
% TR=1.2;
% fout = 0.025;
% for sub = 1:Nsubs
%     dtemp = data{sub};
%     dnew = zeros(size(dtemp));
%     parfor v = 1:size(dtemp,2)
%         dnew(:,v) = filt_notch(spm_vec(dtemp(:,v)),TR,fout,4);
%     end
%     data{sub} = dnew;
% end
% 
% clear dnew dtemp


% 2.1 Normalize and detrend time-series (Takes around 3 seconds per subject).
data = postproc_normalize(data);
inputs.pproc.normalize = 1;


% 2.2 Deconvolve data. This can be done ONLY if data has been normalized.
% define parameters (see and modify within postproc_deconvolve_main.m).

% WARNING: this pproc step is not recommended, but is added as an optional
% step just for cases in which deconvolved data is needed. It is suggested
% to deconvolve first and save the dataset. Also, it takes a long time to
% perform voxel-wise deconvolution (around 15 minutes per subject)

% [dec_out, data] = postproc_deconvolve_main(data);
% inputs.pproc.deconvolve = 1;


% 2.4 Censor frames flagged with high motion (takes 0.5 seconds per subject).
data = postproc_censor_data(data, motion_info);
inputs.pproc.censor = 1;


% Finally assure that data is in single format
for sub = 1:Nsubs
    data{sub} = single(data{sub});
end

disp('postprocessing done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Clustering of data. -> CHOOSE THE DATA OF THE CORRECT POST-PROCESSING STEP TO WORK WITH!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before running the k-means algorithm, define the parameters that will be
% used, and also how to mask data as in Liu et al. [2]

% save the amount of observations for each subject.
for sub = 1:length(data)
    inputs.Nobs(sub) = size(data{sub},1);
end


% parcellate data

% define the seeds.
inputs.seed_dir = '/home/safaai/DanielG/rois_atlas/rois_ludo/';
cd(inputs.seed_dir)

seed_list = dir('*.nii.gz');
% seed_order = [3,6,10,4,5,7,8,1,2,9];
Nrois = length(seed_list);
for s = 1:length(seed_list)
    inputs.seed_names{s} = seed_list(s).name(1:end-4);
    inputs.seed_paths{s} = [pwd '\' seed_list(s).name];
end

%turn data_chosen into cell
tt=1;
for sub = 1:Nsubs
    data_chosen_tmp{sub}(1:inputs.Nobs(sub),:) = data(tt:tt+inputs.Nobs(sub)-1,:);
    tt=tt+inputs.Nobs(sub);
end



%Extract data from seeds and the GS and normalize it. Also compute spectra
%extract seed TS.
inputs.TR = 1; % define sampling rate (s)
flag_roisNAN = zeros(Nsubs,Nrois);
for s = 1:length(inputs.seed_names)
    seed_data{s} = caps_SB_extract_TS(data, inputs.seed_paths{s}, inputs.seed_names{s}, mask_info.mask);    
    seed_data{s}.TS_norm = postproc_normalize(seed_data{s}.TS);    
    seed_data{s}.PSD = caps_analysis_psd(seed_data{s}.TS_norm,inputs.TR);
    tt=1;
    for sub = 1:Nsubs
        Nobs = inputs.Nobs(sub);
        data_chosen_seed1(tt:tt+Nobs-1,s) = seed_data{s}.TS_norm{sub};
        tt = tt+Nobs;
        if sum(seed_data{s}.TS{sub})==0 || sum(isnan(seed_data{s}.TS_norm{sub}))==Nobs
            flag_roisNAN(sub,s)=1;
        end
    end
    
end

% eliminate ROIs with NaN or zeros.
data_chosen_seed=[];
s_on = 0;
for s = 1:Nrois
    if sum(flag_roisNAN(:,s))==0
        s_on=s_on+1;
        data_chosen_seed(:,s_on) = data_chosen_seed1(:,s);
    end
end



end %function