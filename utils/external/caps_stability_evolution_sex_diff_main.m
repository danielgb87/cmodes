


%insert the number of different runs and also the maximum 
N_runs=10;
Max_caps=30;
base_path = '/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/10_runs_corrected';
cd(base_path)
%here the script calculates the stability evolution of caps for TD_M subjects
for i=1:N_runs
    %results_folder = 'MRC_TD_Craddock_950parcel_r1_plus_k2-40';
    results_folder = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/TD_M/TD_M_10_runs_corrected/MRC_TD_M_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
    cd(results_folder)
    caps_parcel_cluster_centroid_stability_evolution(Max_caps);
    cd ..
end

%here the script calculates the stability evolution of caps for TD_F subjects
for i=1:N_runs
    %results_folder = 'MRC_TD_Craddock_950parcel_r1_plus_k2-40';
    results_folder = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/TD_F/TD_F_10_runs_corrected/MRC_TD_F_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
    cd(results_folder)
    caps_parcel_cluster_centroid_stability_evolution(Max_caps);
    cd ..
end

%here the script calculates the stability evolution of caps for AD_M subjects
for i=1:N_runs
    %results_folder = 'MRC_AD_Craddock_950parcel_r1_plus_k2-40';
    results_folder = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/AD_M/AD_M_10_runs_corrected/MRC_AD_M_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
    cd(results_folder)
    caps_parcel_cluster_centroid_stability_evolution(Max_caps);
    cd ..
end

%here the script calculates the stability evolution of caps for AD_F subjects
for i=1:N_runs
    %results_folder = 'MRC_AD_Craddock_950parcel_r1_plus_k2-40';
    results_folder = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/AD_F/AD_F_10_runs_corrected/MRC_AD_F_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
    cd(results_folder)
    caps_parcel_cluster_centroid_stability_evolution(Max_caps);
    cd ..
end

% in this section, the scripts performe the caps stability analysis between
% runs for the same k
% In order to run this part you need to add the following file to your
% path: "caps_parcel_compare_groups_btw_runs_stability.m"
inputs.base_path = base_path;
cd(inputs.base_path)
k=1;
for i=1:N_runs
        runs_TD_M{i} = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/TD_M/TD_M_10_runs_corrected/MRC_TD_M_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
        runs_AD_M{i} = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/AD_M/AD_M_10_runs_corrected/MRC_AD_M_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
end
k_results=cell(Max_caps-1,1);
for t=2:Max_caps
    k_results{k}=caps_parcel_compare_groups_btw_runs_stability(inputs,runs_TD_M,runs_AD_M,N_runs,k);
    k=k+1;
end
save('compare_group_btw_runs_stability_TD_and_AD_M.mat','k_results','-v7.3');

inputs.base_path = base_path;
cd(inputs.base_path)
k=1;
for i=1:N_runs
        runs_TD_F{i} = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/TD_F/TD_F_10_runs_corrected/MRC_TD_F_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
        runs_AD_F{i} = ['/media/DATA/fgatto/MRC_Human_CAPs_analysis/Human_CAPs_Craddock_parcellation/AD_F/AD_F_10_runs_corrected/MRC_AD_F_Craddock_950parcel_r' num2str(i) '_plus_k2-30'];
end
k_results=cell(Max_caps-1,1);
for t=2:Max_caps
    k_results{k}=caps_parcel_compare_groups_btw_runs_stability(inputs,runs_TD_F,runs_AD_F,N_runs,k);
    k=k+1;
end
save('compare_group_btw_runs_stability_TD_and_AD_F.mat','k_results','-v7.3');
