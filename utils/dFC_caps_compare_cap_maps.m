function dFC_caps_compare_cap_maps

clc,clear
%% EDIT HERE
%% 1.a load caps results, inputs, mask, and parcellation info.
% Also define the main directory to save analyses.
k=10;
main = ['/media/DATA/dgutierrez/human_data/caps_human/msc_5sessions_craddock950_ppTingDGB/analysis_parcel/cap_maps_3mm_k' num2str(k)];
cd(main)
ndatasets = 5; % define the amount of datasets to compare
ncaps = 2:20;    % the clusterings you want to compare... make sure all datasets and runs have them.
comp_name = 'caps_parcel_scrub_vs_noScrub';
%Load the results from each datasets for the best run
load('/media/DATA/dgutierrez/human_data/caps_human/msc_5sessions_craddock950_ppTingDGB/analysis_parcel/cap_maps_3mm_k10/cap_results_k10.mat')
caps_concat1 = results_map.map_vw.cap_mean_map;
for nd = 1:ndatasets
    caps_ds1{nd} = results_map.map_vw_ds{nd}.cap_mean_map;
end
clear results_map

load('/media/DATA/dgutierrez/human_data/caps_human/msc_5sessions_craddock950_ppTingDGB/analysis_parcel_scrub/cap_maps_3mm_k10/cap_results_k10.mat')
caps_concat2 = results_map.map_vw.cap_mean_map;
for nd = 1:ndatasets
    caps_ds2{nd} = results_map.map_vw_ds{nd}.cap_mean_map;
end
clear results_map

mkdir(comp_name)
cd(comp_name)

%% compare mean maps (match then compare
for i = 1:k
    for j = 1:k
        Dmat_concat(i,j) = pdist([(spm_vec(caps_concat1(i,:))'); (spm_vec(caps_concat2(j,:))')]);
        for nd = 1:ndatasets
            Dmat_ds{nd}(i,j) = pdist([(spm_vec(caps_ds1{nd}(i,:))'); (spm_vec(caps_ds2{nd}(j,:))')]);
        end
    end
end
[ind_concat, cost_concat] = munkres_HA(Dmat_concat);
for nd = 1:ndatasets
    [ind_ds{nd}, cost_ds{nd}] = munkres_HA(Dmat_ds{nd});
end


% Compute pairwise similarity between matched CAPs from each ds.
for c = 1:k
    cap_sim_concat(c) = corr(spm_vec(caps_concat1(c,:)),spm_vec(caps_concat2(ind_concat(c),:)));
    for c2 = 1:k
        cap_corr_mat_concat(c,c2) = corr(spm_vec(caps_concat1(c,:)),spm_vec(caps_concat2(ind_concat(c2),:)));
    end
end

for nd = 1:ndatasets
    for c = 1:k
        cap_sim_ds{nd}(c) = corr(spm_vec(caps_ds1{nd}(c,:)),spm_vec(caps_ds2{nd}(ind_ds{nd}(c),:)));
        for c2 = 1:k
            cap_corr_mat_ds{nd}(c,c2) = corr(spm_vec(caps_ds1{nd}(c,:)),spm_vec(caps_ds2{nd}(ind_ds{nd}(c2),:)));
        end
    end
end
    
%% 3 plot results
dFC_utils_plot_imagesc_matrix_with_values(cap_corr_mat_concat, 'jet')
print(['cap_corr_concat'], '-dpng','-r600')
close all
for nd = 1:ndatasets
    dFC_utils_plot_imagesc_matrix_with_values(cap_corr_mat_ds{nd}, 'jet')
    print(['cap_corr_ds' num2str(nd)], '-dpng','-r600')
    close all
end

cap_comp.cap_corr_mat_concat = cap_corr_mat_concat;
cap_comp.cap_corr_mat_ds = cap_corr_mat_ds;
cap_comp.cap_sim_concat = cap_sim_concat;
cap_comp.cap_sim_ds = cap_sim_ds;
cap_comp.ind_concat = ind_concat;
cap_comp.ind_ds = ind_ds;
save(comp_name,'cap_comp')




end %function




