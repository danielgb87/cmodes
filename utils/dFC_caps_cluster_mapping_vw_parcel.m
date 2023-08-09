function dFC_caps_cluster_mapping_vw_parcel

% This is a function to cluster and map already processed data. Put it in
% the folder where your save input, data, parcel_data, and mask info is and
% edit the first lines. The script takes the data and parameters, clusters
% data into CAPs at the group level, then makes group and subject level
% CAPs. This is done at the parcel level, or optionally also at the
% voxelwise level.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 26-02-2020
% _________________________________________________________________________


clc,clear
%% EDIT HERE
% 1.a load data, inputs, mask, and parcellation info.
main = '/home/dgutierrez/analyses_daniel/caps_analysis/dbh_kord/gqneg_pre_cno_bs';
cd(main)
load('data_gqneg_pre_cno_bs.mat')
load('data_parcel_gqneg_pre_cno_bs.mat')
load('inputs_gqneg_pre_cno_bs.mat')
load('mask_info_gqneg_pre_cno_bs.mat')
load('subject_list_gqneg_pre_cno_bs.mat')
if isenum(inputs.path_motion_data) == 0
    load('motion_info_gqneg_pre_cno_bs.mat')
end

analysis_id = 'all_subjects';

% 1.b check the motion information, and decide if exclude subjects from
% analysis.
rm_subjects = [];

% 1.c decide which clustering to do and its parameters
do_caps_parcel = 1;
do_caps_vw = 0;


%parcel caps inputs
if do_caps_parcel ==1
    inputs.caps_parcel_nclust = 4:8;    % range of k's (amount of clusters) to do
    inputs.caps_parcel_max_iter = 1000;     % maximum amount of iterations to reach convergence.
    inputs.caps_parcel_Nreps = 10;     % number of initializations of the algorithm.          
    inputs.caps_parcel_distance = 'correlation';   % distance metric to use (see kmeans.m).
    inputs.caps_parcel_online_phase = 'on';    % guarantees a solutio (see kmeans.m)
    inputs.caps_parcel_opts = statset('UseParallel',1);    %CHECK THE AMOUNT OF CORES TO USE SO YOU DON'T BLOW UP THE MACHINE
    inputs.caps_parcel_start = 'plus';      % The initialization condition (see kmeans.m), default use 'plus' or 'centroids'

    % In case you have a preliminary hypothesis of the cluster "centroids", use
    % CAPs, use start = 'centroids' and load a set of hypothesis CAPs (K x N)
    % matrix with the seed-clusters to use. The following part will
    % automatically create seeds and put them in order. Remember that it can
    % only be run for the amount of clusters that the centroids file has.

    if strcmp(inputs.caps_parcel_start, 'centroids') ==1
        if length(inputs.caps_parcel_nclust) > 1
            warning('more than one amount of clusters is done!!!')
        end
        inputs.caps_parcel_centroid_init    = [];
        clusters_init = 'D:\Gutierrez\mice_dFC\results_paper\cap_templates7.mat';
        load(clusters_init)
        inputs.caps_parcel_centroid_init = cap_templates;
    end
    mkdir(['caps_parcel_' analysis_id])
end


% voxelwise CAPs inputs
if do_caps_vw == 1
    inputs.caps_vw_Liu_mask = 1;
    inputs.caps_vw_Liu_mask_top = 0.1;  % highest % of voxels per frame to keep
    inputs.caps_vw_Liu_mask_low = 0.1;  % lowest % of voxels per frame to keep
    inputs.caps_vw_nclust = 4:10;    % range of k's (amount of clusters) to do
    inputs.caps_vw_max_iter = 1000;     % maximum amount of iterations to reach convergence.
    inputs.caps_vw_Nreps = 10;     % number of initializations of the algorithm.          
    inputs.caps_vw_distance = 'correlation';   % distance metric to use (see kmeans.m).
    inputs.caps_vw_online_phase = 'on';    % guarantees a solutio (see kmeans.m)
    inputs.caps_vw_opts = statset('UseParallel',1);    %CHECK THE AMOUNT OF CORES TO USE SO YOU DON'T BLOW UP THE MACHINE
    inputs.caps_vw_start = 'plus';      % The initialization condition (see kmeans.m), default use 'plus' or 'centroids'

    % In case you have a preliminary hypothesis of the cluster "centroids", use
    % CAPs, use start = 'centroids' and load a set of hypothesis CAPs (K x N)
    % matrix with the seed-clusters to use. The following part will
    % automatically create seeds and put them in order. Remember that it can
    % only be run for the amount of clusters that the centroids file has.

    if strcmp(inputs.caps_vw_start, 'centroids') ==1
        if length(inputs.caps_vw_nclust) > 1
            warning('more than one amount of clusters is done!!!')
        end
        inputs.caps_vw_centroid_init    = [];
        clusters_init = 'D:\Gutierrez\mice_dFC\results_paper\cap_templates7.mat';
        load(clusters_init)
        inputs.caps_vw_centroid_init = cap_templates;
    end    
    mkdir(['caps_vw_' analysis_id])
end

%% 2. Organize data, do spatial masking and remove subjects if necessary
nsubs = length(data);
nvox = size(data{1},2);
for sub = 1:nsubs
    nobs(sub) = size(data{sub},1);
    inputs.Nobs(sub) = nobs(sub);
end

if do_caps_vw ==1
    if isempty(rm_subjects) == 0
        data(rm_subjects) = [];
        nsubs = nsubs-length(rm_subjects);
        nobs(rm_subject) = [];
    end
    
    if inputs.caps_vw_Liu_mask == 1
        for sub = 1:nsubs
         data_chosen_vw{sub,1} = dFC_utils_postproc_liu_mask(data{sub}, inputs.caps_vw_Liu_mask_top, inputs.caps_vw_Liu_mask_low);
        end
    end
end

if do_caps_vw ==1
    if isempty(rm_subjects) == 0
        data_chosen_parcel(rm_subjects) = [];
        if do_caps_vw ==0
            nsubs = nsubs-length(rm_subjects);
            nobs(rm_subjects) = [];
        end
    else
        data_chosen_parcel = data_parcel;
    end
end

%% 3. Perform clustering
cd(main)
% parcel level
    data_chosen_parcel = cell2mat(data_chosen_parcel);
    data_chosen_vw = cell2mat(data_chosen_vw);
if do_caps_parcel == 1
    
    nclust = inputs.caps_parcel_nclust;
    iter = inputs.caps_parcel_max_iter;
    dist = inputs.caps_parcel_distance;
    reps = inputs.caps_parcel_Nreps;
    onlinePhase = inputs.caps_parcel_online_phase;
    start = inputs.caps_parcel_start;
    opt = inputs.caps_parcel_opts;
    cd(['caps_parcel_' analysis_id])
    caps_results_parcel = dFC_utils_caps_kmeans(data_chosen_parcel,nsubs,nobs,nclust,dist,iter,reps,onlinePhase,start,opt);
        
        % compute and map CAPs using data (not masked). Do this at the
        % group level and single-subject
        for k = nclust
            caps_results_parcel.(['caps_' num2str(k)]) = dFC_utils_caps_mapping_group_concat(data_chosen_vw,caps_results_parcel.(['caps_' num2str(k)]),k);
            caps_results_parcel.(['caps_' num2str(k)]) = dFC_utils_caps_mapping_ss(data,caps_results_parcel.(['caps_' num2str(k)]),k);
            mkdir (['cap_maps_k' num2str(k)])
            cd (['cap_maps_k' num2str(k)])
            for c = 1:k
                dFC_utils_caps_create_nii(spm_vec(caps_results_parcel.(['caps_' num2str(k)]).cap_mean_map(c,:)), mask_info, (['concat_mean_cap_' num2str(c)]));
                dFC_utils_caps_create_nii(spm_vec(caps_results_parcel.(['caps_' num2str(k)]).cap_T_map(c,:)), mask_info, (['concat_T_cap_' num2str(c)]));
                for sub = 1:nsubs
                    dFC_utils_caps_create_nii(spm_vec(caps_results_parcel.(['caps_' num2str(k)]).cap_mean_map_ss{sub}(c,:)), mask_info, (['ss_mean_cap_' num2str(c) '_s' num2str(sub)]));
                    dFC_utils_caps_create_nii(spm_vec(caps_results_parcel.(['caps_' num2str(k)]).cap_T_map_ss{sub}(c,:)), mask_info, (['ss_T_cap_' num2str(c) '_s' num2str(sub)]));
                end
            end
            cd ..
            
        end
        
    save('caps_results_parcel','caps_results_parcel','-v7.3')
    save('inputs','inputs','-v7.3')
 
end

% voxelwise level
cd(main)
if do_caps_vw == 1
    
    nclust = inputs.caps_vw_nclust;
    iter = inputs.caps_vw_max_iter;
    dist = inputs.caps_vw_distance;
    reps = inputs.caps_vw_Nreps;
    onlinePhase = inputs.caps_vw_online_phase;
    start = inputs.caps_vw_start;
    opt = inputs.caps_vw_opts;
    cd(['caps_vw_' analysis_id])
    caps_results_vw = dFC_utils_caps_kmeans(data_chosen_vw,nsubs,nobs,nclust,dist,iter,reps,onlinePhase,start,opt);
        
        % compute and map CAPs using data (not masked). Do this at the
        % group level and single-subject
        for k = nclust
            caps_results_vw.(['caps_' num2str(k)]) = dFC_utils_caps_mapping_group_concat(data_chosen_vw,caps_results_vw.(['caps_' num2str(k)]),k);
            caps_results_vw.(['caps_' num2str(k)]) = dFC_utils_caps_mapping_ss(data,caps_results_vw.(['caps_' num2str(k)]),k);
            mkdir (['cap_maps_k' num2str(k)])
            cd (['cap_maps_k' num2str(k)])
            for c = 1:k
                dFC_utils_caps_create_nii(spm_vec(caps_results_vw.(['caps_' num2str(k)]).cap_mean_map(c,:)), mask_info, (['concat_mean_cap_' num2str(c)]));
                dFC_utils_caps_create_nii(spm_vec(caps_results_parcel.(['caps_' num2str(k)]).cap_T_map(c,:)), mask_info, (['concat_T_cap_' num2str(c)]));
                for sub = 1:nsubs
                    dFC_utils_caps_create_nii(spm_vec(caps_results_vw.(['caps_' num2str(k)]).cap_mean_map_ss{sub}(c,:)), mask_info, (['ss_mean_cap_' num2str(c) '_s' num2str(sub)]));
                    dFC_utils_caps_create_nii(spm_vec(caps_results_vw.(['caps_' num2str(k)]).cap_T_map_ss{sub}(c,:)), mask_info, (['ss_T_cap_' num2str(c) '_s' num2str(sub)]));
                end
            end
            cd ..
            
        end
        
    save('caps_results_vw','caps_results_vw','-v7.3')
    save('inputs','inputs','-v7.3')
    
        
end


disp('done with clustering')

% end main function.
end



