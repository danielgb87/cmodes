function caps_05_dFC_caps_cluster_map

% This is a function to map already clustered frames. you can do it using
% results from any analysis (parcellated, vw or vw-downsampled) but map the
% averaged frames of another, for example use results from parcellated CAPs
% to map their centroids by averaging voxelwise frames. To do this, you
% have to specify if your clusters are in the same source as the data you
% are loading.

%Here you can also decide if map the CAPs also at the subject level. !!!
%make sure that you are mapping from the best run for each K, or else run
%this script separately by modifying the fulder of CAPs results with the
%best run for that specific K.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 01-27-2021
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% path to dFC_CAPs_scripts and external toolboxes
addpath(genpath('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127'))

% 1.a load data, inputs, mask, and parcellation info.
main = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221';
addpath(genpath(main))

main_caps = [main '/caps_vw_awk_ane_plus_r2']; %folder of caps_results
main_analysis = [main '/analysis_vw_plus']; %directory where to save

% load data to map, its mask and subject list, remember to change the name
% of data
load('data_vw_awk_lowres.mat')
data = data_vw;
clear data_vw
load('mask_info_awk_lowres.mat')

% 1.b load caps results from the run of your choice, remember to change the
% name, and select configurations (ncaps) that were computed in your results
cd(main_caps)
ncaps = [6]; %this can be a vector with diffrent k's. One folder will be create
load('inputs.mat')
load([inputs.analysis_name])
ds_name=inputs.ds_name;

%load the results from K of interest, then clear the results structure
for k = 1:length(ncaps)
    results{k} = caps_results.(['caps_' num2str(ncaps(k))]); % CHANGE NAME HERE!!!
end
clear caps_results;

% 1.c define if the data is the same one used for the clustering. Elsewise,
% the algorithm will recompute the CAP maps for the chosen data.
% caps_from_source_data=1 means themaps will be done from the same space as
% the clustering results, =0 means that they have to be recomputed based on
% the clustering results.
caps_from_source_data =1;
data_source = 'lowres';

% decide which caps to map and if do you want to map single-subject CAPs
% (yes=1, no=0)
do_ss = 0;

% remove subjects and organize data
rm_subjects = inputs.rm_subjects;
ss=0;
for nd = 1:length(data)
    nsubs(nd) = length(data{nd});
    if isempty(rm_subjects{nd}) == 0
        data{nd}(rm_subjects{nd}) = [];
        nsubs(nd)=nsubs(nd)-length(rm_subjects{nd});
    end
    data_mat{nd} = cell2mat(data{nd});
end
data_mat = cell2mat(data_mat');

for nd = 1:length(data)
    for sub = 1:length(data{nd})
        ss=ss+1;
        data_concat_cell{ss} = data{nd}{sub};
        sub_list_ds{nd}(sub) = ss;
        nobs_list_ds{nd}(sub) = size(data{nd}{sub},1);
    end
end
ndatasets = length(data);
clear data
nsubs_concat=ss;

%% 2. map CAPs        
% group level and single-subject if do_ss=1
cd(main_analysis)
for k = 1:length(ncaps)
    mkdir (['cap_maps_' data_source '_k' num2str(ncaps(k))])
    cd (['cap_maps_' data_source '_k' num2str(ncaps(k))])
    if caps_from_source_data==0
        results_new{k}.frame_index = results{k}.frame_index;
        results_new{k} = dFC_utils_caps_mapping_group_concat(data_mat,results_new{k},ncaps(k));
        % now do caps for each dataset
        nn=1;
        for nd=1:ndatasets
            results_new_ds{nd,k}.frame_index = results{k}.frame_index(nn:nn+sum(nobs_list_ds{nd}-1));
            results_new_ds{nd,k} = dFC_utils_caps_mapping_group_concat(data_mat(nn:nn+sum(nobs_list_ds{nd})-1,:),results_new_ds{nd,k},ncaps(k));
            nn = nn+sum(nobs_list_ds{nd});
        end
        for c = 1:ncaps(k)
            dFC_utils_caps_create_nii(spm_vec(results_new{k}.cap_mean_map(c,:)), mask_info, (['concat_mean_cap_' num2str(c)]));
            dFC_utils_caps_create_nii(spm_vec(results_new{k}.cap_T_map(c,:)), mask_info, (['concat_T_cap_' num2str(c)]));
            for nd = 1:ndatasets
                dFC_utils_caps_create_nii(spm_vec(results_new_ds{nd,k}.cap_mean_map(c,:)), mask_info, ([ds_name{nd} '_mean_cap_' num2str(c)]));
                dFC_utils_caps_create_nii(spm_vec(results_new_ds{nd,k}.cap_T_map(c,:)), mask_info, ([ds_name{nd} '_T_cap_' num2str(c)]));
            end
            if do_ss ==1
                results_new{k}.frame_ind_sub = results{k}.frame_ind_sub;
                results_new{k} = dFC_utils_caps_mapping_ss(data_concat_cell,results_new{k},ncaps(k));

                for nd = 1:ndatasets
                    for sub = 1:nsubs(nd)
                        dFC_utils_caps_create_nii(spm_vec(results_new{k}.cap_mean_map_ss{sub_list_ds{nd}(sub)}(c,:)), mask_info, (['ss_mean_cap_' num2str(c) '_' ds_name{nd} '_s' num2str(sub)]));
                        dFC_utils_caps_create_nii(spm_vec(results_new{k}.cap_T_map_ss{sub_list_ds{nd}(sub)}(c,:)), mask_info, (['ss_T_cap_' num2str(c) '_' ds_name{nd} '_s' num2str(sub)]));
                    end
                end
            end
        end
        results_map = results{k};
        results_map.map_vw = results_new{k};
        for nd = 1:ndatasets
            results_map.map_vw_ds{nd} = results_new_ds{nd,k};
        end
        save(['cap_results_k' num2str(ncaps(k))] ,'results_map','-v7.3')
        save('inputs' ,'inputs','-v7.3')
        clear results_map
    else
        %data and CAPs come from the same source, no need to recompute maps
        %  do caps for each dataset
        nn=1;
        for nd=1:ndatasets
            results_ds{nd,k}.frame_index = results{k}.frame_index(nn:nn+sum(nobs_list_ds{nd})-1);
            results_ds{nd,k}= dFC_utils_caps_mapping_group_concat(data_mat(nn:nn+sum(nobs_list_ds{nd})-1,:),results_ds{nd,k},ncaps(k));
            nn = nn+sum(nobs_list_ds{nd});
        end
        for c = 1:ncaps(k)
            dFC_utils_caps_create_nii(spm_vec(results{k}.cap_mean_map(c,:)), mask_info, (['concat_mean_cap_' num2str(c)]));
            dFC_utils_caps_create_nii(spm_vec(results{k}.cap_T_map(c,:)), mask_info, (['concat_T_cap_' num2str(c)]));
            for nd = 1:ndatasets
                dFC_utils_caps_create_nii(spm_vec(results_ds{nd,k}.cap_mean_map(c,:)), mask_info, ([ds_name{nd} '_mean_cap_' num2str(c)]));
                dFC_utils_caps_create_nii(spm_vec(results_ds{nd,k}.cap_T_map(c,:)), mask_info, ([ds_name{nd} '_T_cap_' num2str(c)]));
            end
            if do_ss ==1
                for nd = 1:ndatasets
                    for sub = 1:nsubs(nd)
                        dFC_utils_caps_create_nii(spm_vec(results{k}.cap_mean_map_ss{sub_list_ds{nd}(sub)}(c,:)), mask_info, (['ss_mean_cap_' num2str(c) '_' ds_name{nd} '_s' num2str(sub)]));
                        dFC_utils_caps_create_nii(spm_vec(results{k}.cap_T_map_ss{sub_list_ds{nd}(sub)}(c,:)), mask_info, (['ss_T_cap_' num2str(c) '_' ds_name{nd} '_s' num2str(sub)]));
                    end
                end
            end
        end
        results_map = results{k};
        results_map.map_vw = results{k};
        for nd = 1:ndatasets
            results_map.map_vw_ds{nd} = results_ds{nd,k};
        end
        save(['cap_results_k' num2str(ncaps(k))] ,'results_map','-v7.3')
        save('inputs' ,'inputs','-v7.3')
        clear results_map
        
    end
    cd ..
end






end % function