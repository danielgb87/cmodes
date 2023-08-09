function dFC_FC_seed_regress_cap

% This function takes two rois, their timecourses, and FC data, and  using
% the CAPs information, plots the timeocurse at the occurrence of each CAP.
% It also computes the dFC at each CAP, and recomputes tghe FC after
% regressing each CAP timecourse.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 7-03-2020
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% 1.a load data.
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC_parcel/SBA/FC_seeds';
main_data = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC_parcel/SBA';
cd(main_data)
load('data_parcel_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('inputs_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('subject_list_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
ndatasets = 3;
for nd = 1:ndatasets
    nsubs(nd) = length(data_parcel{nd});
end
nrois = inputs.Nrois;

%% 2. Compute FC for each dataset
cd(main)






end % function