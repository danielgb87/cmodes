function dFC_SBA_FC_mCAPs


% This function takes extracted data and a set of seeds,a nd computes for
% each group the FC and mCAP maps for each seed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Daniel Gutierrez-Barragan 2020, V1 10-04-20
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% 1.a load data.
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC/SBA_vw_ds/';
main_data = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/';
seed_path = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC/SBA_vw_ds/seeds';
cd(main_data)
load('data_vw_scrub_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
scrub = 1;
data = data_vw_scrub;
clear data_vw_scrub
load('inputs_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('subject_list_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('mask_info_mac_concat_2mm_woVent_WM_Cereb_BS')
ndatasets = 3;
for nd = 1:ndatasets
    nsubs(nd) = length(data{nd});
end
cd(seed_path)
seed_list = dir('*.nii.gz');
%extract TS of seeds
for nd = 1:ndatasets
    for s = 1:length(seed_list)
        [seed_data{nd,s}] = dFC_utils_SBA_ts_extract(data{nd}, seed_list(s).name, mask_info.mask);
    end
end
nseeds = length(seed_list);

%% 2. compute FC, mCAPs, and their comparison
cd(main)
% Compute the connectivity maps
mkdir('subject_maps')

for s = 1:nseeds
    for nd = 1:ndatasets
        [seed_data{nd,s}.CMaps, seed_data{nd,s}.CMaps_mean, seed_data{nd,s}.CMaps_T] = dFC_utils_SBA_FCcomp(data{nd},seed_data{nd,s}.TS);
        [seed_data{nd,s}.FC_Params] = dFC_utils_SBA_FCparams(seed_data{nd,s}.TS, data{nd}, seed_data{nd,s}.CMaps); 
        seed_data{nd,s}.seed_thr = 15; % percent of high BOLD frames to use.
        [seed_data{nd,s}.selected_frames, seed_data{nd,s}.frame_list] = dFC_utils_SBA_seed_select_frames(seed_data{nd,s}.TS, seed_data{nd,s}.seed_thr,seed_data{nd,s}.FC_Params);
        
        for sub = 1:nsubs(nd)
            nobs(sub) = size(data{nd}{sub},1);
            data_chosen_seed{nd,s}{sub} = data{nd}{sub}(seed_data{nd,s}.frame_list{sub},:);
            seed_data{nd,s}.mCAP_ss{sub} = mean(data_chosen_seed{nd,s}{sub},1);
        end
        dtemp=cell2mat(data_chosen_seed{nd,s}');
        seed_data{nd,s}.mCAP_concat = mean(dtemp,1);
        for v = 1:size(dtemp,2)
            [~,~,~,t] = ttest(spm_vec(dtemp(:,v)),0,'tail','both');
            seed_data{nd,s}.mCAP_concat_T(v) = t.tstat;     
        end
        cd('subject_maps')
        for sub = 1:nsubs(nd)
            seed_data{nd,s}.nobs_seed(sub) = size(data_chosen_seed{nd,s}{sub},1);
            dFC_utils_caps_create_nii(seed_data{nd,s}.CMaps{sub}, mask_info, ['FC_s' num2str(sub) '_ds' num2str(nd) '_' seed_list(s).name(1:end-7)])
            dFC_utils_caps_create_nii(seed_data{nd,s}.mCAP_ss{sub}, mask_info, ['mCAP_s' num2str(sub) '_ds' num2str(nd) '_' seed_list(s).name(1:end-7)])
        end
        
       cd(main)

    end
    
end

%plot the FC vs mCAP similarity curves
for s = 1:nseeds
    for nd = 1:ndatasets
        dFC_utils_SBA_plot_seedFC2mCAP_corr(['FC2cap_corr_ds' num2str(nd) '_' seed_list(s).name(1:end-7)], seed_data{nd,s}.FC_Params);
        close all          
    end
end

%save seed_data
save('seed_data','seed_data')


% end main function.
end



