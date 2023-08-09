function dFC_FC_global_connectivity


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 26-02-2020
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% 1.a load data.
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC_parcel/';
main_data = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled';
cd(main_data)
load('data_vw_scrub_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('inputs_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('subject_list_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('mask_info_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
data = data_vw_scrub; clear data_vw_scrub
ndatasets = 3;
for nd = 1:ndatasets
    nsubs(nd) = length(data{nd});
end
nrois = inputs.Nrois;

%% 2 compute global connectivity for all subjects, including negative values of correaltion
cd(main)
mkdir('global_FC')
cd('global_FC')
for nd = 1:ndatasets
    for sub = 1:nsubs(nd)
        [nobs, nvox] = size(data{nd}{sub});        
        ds = data{nd}{sub};
        fc = atanh(corr(ds));
        fc(isnan(fc))=0;
        fc(isinf(fc))=0;
        fc_abs = abs(fc);
        gFCz{nd}(sub,:) = mean(fc,1);
        gFCz_abs{nd}(sub,:) = mean(fc_abs,1);
    end
    gFC_mean{nd} = mean(gFCz{nd},1);
    gFC_abs_mean{nd} = mean(gFCz_abs{nd},1);
    for v = 1:nvox
        [gFC_h{nd}(v),gFC_p{nd}(v),~,t] = ttest(spm_vec(gFCz{nd}(:,v)),0,'tail','both');
        gFC_T{nd}(v) = t.tstat;
        [gFC_abs_h{nd}(v),gFC_abs_p{nd}(v),~,t] = ttest(spm_vec(gFCz_abs{nd}(:,v)),0,'tail','right');
        gFC_abs_T{nd}(v) = t.tstat;
    end
    dFC_utils_caps_create_nii(spm_vec(tanh(gFC_mean{nd})), mask_info, ['gFC_mean_ds' num2str(nd)])
    dFC_utils_caps_create_nii(spm_vec(tanh(gFC_abs_mean{nd})), mask_info, ['gFC_abs_mean_ds' num2str(nd)])
    dFC_utils_caps_create_nii(spm_vec(gFC_T{nd}), mask_info, ['gFC_T_ds' num2str(nd)])
    dFC_utils_caps_create_nii(spm_vec(gFC_abs_T{nd}), mask_info, ['gFC_abs_T_ds' num2str(nd)])
end

%% 3. compare maps between datasets
for nd1= 1:ndatasets-1
    for nd2=nd1+1:ndatasets
        for v =1:nvox
            
            [gFC_diff_h{nd1,nd2}(v), gFC_diff_p{nd1,nd2}(v),~,t] = ttest2(spm_vec(gFCz{nd1}(:,v)),spm_vec(gFCz{nd2}(:,v)),'tail','both');
            gFC_diff_T{nd1,nd2}(v) = t.tstat;
            [gFC_abs_diff_h{nd1,nd2}(v), gFC_abs_diff_p{nd1,nd2}(v),~,t] = ttest2(spm_vec(gFCz_abs{nd1}(:,v)),spm_vec(gFCz_abs{nd2}(:,v)),'tail','both');
            gFC_abs_diff_T{nd1,nd2}(v) = t.tstat;
        end
        
        gFC_diff_T{nd1,nd2}(isinf(gFC_diff_T{nd1,nd2})|isnan(gFC_diff_T{nd1,nd2})) = 0;
        gFC_diff_p{nd1,nd2}(isinf(gFC_diff_p{nd1,nd2})|isnan(gFC_diff_p{nd1,nd2})) = 0;
        
        gFC_abs_diff_T{nd1,nd2}(isinf(gFC_abs_diff_T{nd1,nd2})|isnan(gFC_abs_diff_T{nd1,nd2})) = 0;
        gFC_abs_diff_p{nd1,nd2}(isinf(gFC_abs_diff_p{nd1,nd2})|isnan(gFC_abs_diff_p{nd1,nd2})) = 0;
       
        dFC_utils_caps_create_nii(spm_vec(gFC_diff_T{nd1,nd2}), mask_info, ['gFC_diff_T_ds' num2str(nd1) '-ds' num2str(nd2)])
        dFC_utils_caps_create_nii(spm_vec(gFC_abs_diff_T{nd1,nd2}), mask_info, ['gFC_abs_diff_T_ds' num2str(nd1) '-ds' num2str(nd2)])
    end
end
    











end %function