function caps_02_dFC_caps_cluster_concat

% This function clusters processed and prepared data, concatenating all
% groups. Inputs are the folder where your save input, data, parcel_data,
% and mask info is and edit the first lines. The script takes the data and
% parameters, clusters data into CAPs at the group level, then makes group
% and subject level CAPs.

%This is done at the parcel level level, or optionally also at the
%voxelwise level, depending on input data. The script can also optionally
%use centroid priors from previous studies, or establish the cluster
%centroids for different amounts of k. 

%Be aware that the output is exclusively dependent on the dimensionality of
%the data. Parcellated data will yield results with the amount of parcels,
%not voxelwise. In the next step of mapping you can map voxelwise. Also, if
%you have previously defined centroids, make sure they are on the same
%template of your dta.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V2 27-01-2021
% _________________________________________________________________________

clc,clear

%% 1. EDIT HERE
% path to dFC_CAPs_scripts and external toolboxes
addpath(genpath('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127'))

% 1.a load data, inputs, mask, and parcellation info.
main = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221'; % main folder where data, inputs, mask_info is
addpath(genpath(main))
cd(main)
% 1.b load data, inputs, mask_info, subject_list, and motion info from previous
% step
load('data_vw_awk_lowres.mat')
load('inputs_awk_lowres.mat')
load('mask_info_awk_lowres.mat')
load('subject_list_awk_lowres.mat')
if isenum(inputs.path_motion_data) == 0
    load('motion_info_0750_awk_lowres.mat')
end
data = data_vw; % EDIT THIS WITH THE NAME OF THE DATA STRUCTURE
clear data_vw
inputs.ds_name=inputs.ds_name;

% if scrubbing was done, do it load the scrubbed .mat file. YOU CAN CHOOSE
% IF EITHER USE SCRUBBED OR UNSCRUBBED DATA HERE.

inputs.pproc_scrub = 0;

if inputs.pproc_scrub ==1
    clear data_vw  % EDIT THIS WITH THE NAME OF THE DATA STRUCTURE
    load('data_vw_scrub_awk_lowres.mat')
    data = data_vw_scrub;
    clear data_vw_scrub
end

% 1.c define the amount of datasets (groups/sessions/experimental phases)
% and the ID of your analysis. It is recommended to name the runs also, and
% include in the ID if the analysis is with voxelwise (vw) or parcellated
% data; if you scrubbed; and are using centroids

ndatasets = length(data);
inputs.analysis_name = 'caps_vw_awk_ane_plus_r2';

% 1.d check the motion information, and decide if exclude subjects from
% analysis. Do it for each dataset
for nd = 1:ndatasets
    rm_subjects{nd} = [];
end
inputs.rm_subjects = rm_subjects;

% 1.e decide which clustering to do and its parameters

% voxelwise CAPs inputs... DO NOT USE LIU_MASK IF NOT USING VOXELWISE DATA
inputs.caps_Liu_mask = 1;  % mask frames like Liu et al 2013
inputs.caps_Liu_mask_top = 0.1;  % highest % of voxels per frame to keep
inputs.caps_Liu_mask_low = 0.1;  % lowest % of voxels per frame to keep

% if you are doing CAPs with previous centroids, inputs.caps_Nreps=1, no
% need for more. If you are trying diffrent k-numbers, then do at least 5
% different replications. Iterations for rsfMRI data between 500-1000.
inputs.caps_nclust = [4:6];    % range of k's (amount of clusters) to do
inputs.caps_max_iter = 200;     % maximum amount of iterations to reach convergence.
inputs.caps_Nreps = 1;     % number of initializations of the algorithm.

% decide on the distance metric; if to perform online phase (recommended,
% see kmeans.m); use parallel computing for different replications; and if
% centroids or plus algorithm is used
inputs.caps_distance = 'correlation';   % distance metric to use (see kmeans.m).
inputs.caps_online_phase = 'on';    % guarantees a solutio (see kmeans.m)
inputs.caps_opts = statset('UseParallel',1);    %CHECK THE AMOUNT OF CORES TO USE SO YOU DON'T BLOW UP THE MACHINE
inputs.caps_start = 'plus';      % The initialization condition (see kmeans.m), default use 'plus' or 'centroids'

% In case you have a preliminary hypothesis of the cluster "centroids", use
% CAPs, use inputs.caps_vw_start = 'centroids'above and load a set of CAP
% maps, in a folder where they are located and in the same resolution as
% the input data. The following part will automatically create seeds and
% put them in order. Remember that if centroids are used, it can only be
% run for the amount of clusters that the centroids file has. THE INPUT CAP
% CENTROIDS MUST BE A MERGED FILE WITH A CENTROID IN EACH 'TIMEPOINT'!!!

if strcmp(inputs.caps_start, 'centroids') ==1
    if length(inputs.caps_nclust) > 1
        warning('more than one amount of clusters is done!!!')
    end
    % path to the merged prior CAP files
    inputs.caps_centroid_init    = [];
    [cc, ~, ~] = dFC_utils_prepare_dataset_matrix('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127/cap_centroids_CBpaper_k6/'...
        , 'merged_caps_mean_maps.nii.gz', inputs.path_mask);
    
    inputs.caps_centroid_init = cc{1};

else
    inputs.caps_centroid_init = 'plus';
end

%% 2. Organize data, do spatial masking and remove subjects if necessary
for nd = 1:ndatasets
    nsubs(nd) = length(data{nd});
    nvox = size(data{nd}{1},2);
    for sub = 1:nsubs(nd)
        nobs{nd}(sub) = size(data{nd}{sub},1);
        inputs.nobs{nd}(sub) = nobs{nd}(sub);
    end
end

%mask and remove subjects from voxelwise data
for nd = 1:ndatasets
    if isempty(rm_subjects{nd}) == 0
        data{nd}(rm_subjects{nd}) = [];
        nsubs(nd) = nsubs(nd)-length(rm_subjects{nd});
        nobs{nd}(rm_subjects{nd}) = [];
    end
    if inputs.caps_Liu_mask == 1
        for sub = 1:nsubs(nd)
            data_chosen{nd}{sub,1} = dFC_utils_postproc_liu_mask(data{nd}{sub}, inputs.caps_Liu_mask_top, inputs.caps_Liu_mask_low);
        end
        data_chosen{nd} = cell2mat(data_chosen{nd});
        
    else
        data_chosen{nd} = data{nd};
        data_chosen{nd} = cell2mat(data_chosen{nd});
    end
end

data_concat = cell2mat(data_chosen');

% concatenate data
ss=0;
for nd = 1:ndatasets
    for sub = 1:nsubs(nd)
        ss=ss+1;
        data_concat_cell{ss} = data{nd}{sub};
    end
end
nsubs_concat = sum(nsubs);
nobs_concat = cell2mat(nobs);

%% 3. Perform clusterinG
cd(main)
%make analysis directory
mkdir([inputs.analysis_name])

nclust = inputs.caps_nclust;
iter = inputs.caps_max_iter;
dist = inputs.caps_distance;
reps = inputs.caps_Nreps;
onlinePhase = inputs.caps_online_phase;
start = inputs.caps_centroid_init;
opt = inputs.caps_opts;
cd([inputs.analysis_name])
caps_results = dFC_utils_caps_kmeans(data_concat,nsubs_concat,nobs_concat,nclust,dist,iter,reps,onlinePhase,start,opt);

% compute and map CAPs using data (not Liu-masked). Do this at the group level
% and single-subject
for k = nclust
    caps_results.(['caps_' num2str(k)]) = dFC_utils_caps_mapping_group_concat(data_concat,caps_results.(['caps_' num2str(k)]),k);
    caps_results.(['caps_' num2str(k)]) = dFC_utils_caps_mapping_ss(data_concat_cell,caps_results.(['caps_' num2str(k)]),k);    
end

save([inputs.analysis_name],'caps_results','-v7.3')
save('inputs','inputs','-v7.3')

disp('done with clustering')

% end main function.
end



