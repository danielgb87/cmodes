function dFC_caps_FCpair_BOLD_at_cap


% This function takes a pair of seeds, extract the data and: 1= computes
% their FC in each group and their differences without and after regressing
% each CAP from a set of results. It then constructs the mean BOLD kernel
% at that CAP for each region.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Daniel Gutierrez-Barragan 2020, V1 10-04-20
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% 1.a load data.
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC/SBA_vw_ds/SBA_at_cap/23_24_nba_stp_v1/';
main_data = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/';
seed_path = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC/SBA_vw_ds/SBA_at_cap/23_24_nba_stp_v1/';
cd(main_data)
scrub =1;
rm_subjects =[];
k=8;
% load the data to compute VE from, and rename it
load('data_vw_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
data = data_vw;
clear data_vw;
if scrub==1
    clear data
    load('data_vw_scrub_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
    data = data_vw_scrub;
    clear data_vw_scrub;
end
load('motion_info_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('mask_info_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('subject_list_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('mask_info_mac_concat_2mm_woVent_WM_Cereb_BS')
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_ds/cap_maps_2mm_k8/cap_results_k8')
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_ds/cap_maps_2mm_k8/inputs')
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_ds/cap_maps_2mm_k8/caps_dynamics_k8/caps_dyn')
TR=inputs.analysis_TR;
ndatasets = 3;
for nd = 1:ndatasets
    nsubs(nd) = length(data{nd});
end

% in case subjects were removed
rm_subjects = inputs.rm_subjects;
ss=0;
for nd = 1:ndatasets
    nsubs(nd) = length(data{nd});
    if isempty(rm_subjects{nd}) == 0
        data{nd}(rm_subjects{nd}) = [];
        nsubs(nd)=nsubs(nd)-length(rm_subjects{nd});
        motion_info{nd}.FD_flag{sub} = [];
    end
    data_mat{nd} = cell2mat(data{nd});    
    for sub = 1:nsubs(nd)
        ss=ss+1;
        data_concat_cell{ss} = data{nd}{sub};
        sub_list_ds{nd}(sub) = ss;
    end
end
%correct the FD flags if data was scrubbed
if scrub ==1
   for nd = 1:ndatasets
       for sub = 1:nsubs(nd)
        motion_info{nd}.FD_flag{sub} = zeros(size(data{nd}{sub},1),1);
       end
   end
end

%extract TS of seeds
cd(seed_path)
seed_list = dir('*.nii.gz');
for nd = 1:ndatasets
    for s = 1:length(seed_list)
        [seed_data{nd,s}] = dFC_utils_SBA_ts_extract(data{nd}, seed_list(s).name, mask_info.mask);        
    end
    for sub=1:nsubs(nd)
            nobs{nd}(sub) = length(seed_data{nd,1}.TS{sub});
    end
end
nseeds = length(seed_list);
cd(main)
%% 2. compute and sort peaks of CAP activity
%define the threshold for peaks and the range to compute CAP evo in TRs
cfc_thr = 1;
range = 15;
for nd = 1:ndatasets
    [peak_ind{nd}, peak_amp{nd}] = dFC_utils_caps_cfc_peaks(caps_dyn.cfc.cfc{nd}.cfc_norm, caps_dyn.cfc.cfc{nd}.cfc, motion_info{nd}, cfc_thr);   
end
% take the highest 10% (or change the percentile) of the peaks for each CAP
% in each subject.
top_percent = 0.2;
for nd = 1:ndatasets
    for c = 1:k
        for sub = 1:nsubs(nd)
            nsamples{nd}(c,sub) = ceil(length(peak_ind{nd}{sub,c})*top_percent);
            [A{nd}{sub,c},I{nd}{sub,c}] = sort(peak_amp{nd}{sub,c},'descend');
            peak_ind_ord{nd}{sub,c} = peak_ind{nd}{sub,c}(I{nd}{sub,c}(1:nsamples{nd}(c,sub)));
        end
    end
end

% take he preceding and subsequent "range" timepoints (or change the range) of
% each CAP-peak event, and compute the each ROi's kernel for that range.
% Average them for each subject, then do the group level mean and
% SEM. THE CODE ENSURES THAT THE PEAK IS NOT IN THE BEGGINING OR END OF THE
% TIMECOURSE, DEPENDING ON THE LAGS.
lags = -range:range;

for nd = 1:ndatasets
    for r = 1:nseeds
        for c = 1:k
            
            pk2=0;
            for sub = 1:nsubs(nd)
                [~,~,ts_regc{nd,r}{c}{sub}] = regress(spm_vec(seed_data{nd,r}.TS{sub}),[spm_vec(ones(length(caps_dyn.cfc.cfc{nd}.cfc_norm{sub,c}),1)) spm_vec(caps_dyn.cfc.cfc{nd}.cfc_norm{sub,c})]);

                for pk = 1:length(peak_ind_ord{nd}{sub,c})
                    if peak_ind_ord{nd}{sub,c}(pk)+lags(end) < nobs{nd}(sub) && peak_ind_ord{nd}{sub,c}(pk)+lags(1) > 0
                        pk2 = pk2+1;
                        roi_kernel_at_cap{nd,r}{c}(pk2,:) = seed_data{nd,r}.TS{sub}(peak_ind_ord{nd}{sub,c}(pk)+lags);
                    end
                end
%                 roi_kernel_at_cap_subm{nd}{r,c}(sub,:) = mean(roi_kernel_at_cap{nd}{r,c}{sub},1);
            end                    
            mean_group_roi_kernel{nd,r}(c,:) = mean(roi_kernel_at_cap{nd,r}{c},1);       
            std_group_roi_kernel{nd,r}(c,:) = std(roi_kernel_at_cap{nd,r}{c},1);                       
        end
    end
end

%% 3.
% take each interval of peak-events and compute the correlations between
% the SEEDs for each event and average them, first at the subject level,
% then compuyte the mean and SEM at the group level.

for nd = 1:ndatasets
    for c = 1:k
        pk2=0;
        for sub = 1:nsubs(nd)
            for pk = 1:length(peak_ind_ord{nd}{sub,c})
                if peak_ind_ord{nd}{sub,c}(pk)+lags(end) < nobs{nd}(sub) && peak_ind_ord{nd}{sub,c}(pk)+lags(1) > 0
                    pk2=pk2+1;
                    for r1 = 1:nseeds
                        for r2 = 1:nseeds
                            FC_at_cap_peak{nd,c}(pk2,r1,r2) = corr(spm_vec(roi_kernel_at_cap{nd,r1}{c}(pk2,:)),spm_vec(roi_kernel_at_cap{nd,r2}{c}(pk2,:)));
                        end
                    end
                end
            end
            
            for r1 = 1:nseeds
                for r2=1:nseeds
                    FC{nd}(sub,r1,r2) = corr(spm_vec(seed_data{nd,r1}.TS{sub}), spm_vec(seed_data{nd,r2}.TS{sub}));
                    FC_regc{nd,c}(sub,r1,r2) = corr(spm_vec(ts_regc{nd,r1}{c}{sub}), spm_vec(ts_regc{nd,r2}{c}{sub}));
                end
            end
        end

        for r1 = 1:nseeds
            for r2 = 1:nseeds
                FC_group_mean_at_cap{nd}(r1,r2,c) = mean(FC_at_cap_peak{nd,c}(:,r1,r2));
                FC_group_std_at_cap{nd}(r1,r2,c) = std(FC_at_cap_peak{nd,c}(:,r1,r2));
                FC_group_mean_regc{nd}(r1,r2,c) = mean(FC_regc{nd,c}(:,r1,r2));
                FC_group_std_regc{nd}(r1,r2,c) = std(FC_regc{nd,c}(:,r1,r2));
            end
        end
    end
    for r1 = 1:nseeds
            for r2 = 1:nseeds
                FC_mean{nd}(r1,r2) = mean(FC{nd}(:,r1,r2));
                FC_std{nd}(r1,r2) = std(FC{nd}(:,r1,r2));
            end
    end
        
end
      
% perform stats on local FC
for nd1=1:ndatasets
    for nd2=1:ndatasets
        for c = 1:k
            for r1 = 1:nseeds
                for r2 = 1:nseeds
                    [FC_diff_at_cap_h{nd1,nd2}(r1,r2,c), FC_diff_at_cap_p{nd1,nd2}(r1,r2,c),~ , ~] = ttest2(spm_vec(FC_at_cap_peak{nd1,c}(:,r1,r2)), spm_vec(FC_at_cap_peak{nd2,c}(:,r1,r2)),'tail','both');
                    [FC_diff_regc_h{nd1,nd2}(r1,r2,c), FC_diff_regc_p{nd1,nd2}(r1,r2,c),~ , ~] = ttest2(spm_vec(FC_regc{nd1,c}(:,r1,r2)), spm_vec(FC_regc{nd2,c}(:,r1,r2)),'tail','both');
                    [FC_diff_h{nd1,nd2}(r1,r2), FC_diff_p{nd1,nd2}(r1,r2),~ , ~] = ttest2(spm_vec(FC{nd1}(:,r1,r2)), spm_vec(FC{nd2}(:,r1,r2)),'tail','both');
                end
            end
        end
    end
end


%% 4. plot the ROI kernels and the bargraphs of full and local FC

% plot each activation function for each roi at each CAP
for nd = 1:ndatasets
    for c = 1:k
        for r = 1:nseeds
            figure
            shadedErrorBar(TR*lags,mean_group_roi_kernel{nd,r}(c,:),std_group_roi_kernel{nd,r}(c,:)/sqrt(nsubs(nd)),'k',0.5)
            ax=gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            box off
            set(gca,'color','none') 
            print([seed_list(r).name(1:end-7) '_ds' num2str(nd) '_BOLD_at_cap_' num2str(c)], '-dpng','-r800')       
            close all
        end
    end
end

% plot FC for the whole ts and at each CAP
for nd = 1:ndatasets
    dslabel{nd} = ['ds' num2str(nd)];
end
for s1 = 1:nseeds-1
    s1name = seed_list(s1).name(1:end-7);
    for s2 = s1+1:nseeds
        s2name = seed_list(s2).name(1:end-7);
        
        for nd = 1:ndatasets
            y(1,nd) = FC_mean{nd}(s1,s2);
            e(1,nd) = FC_std{nd}(s1,s2)/sqrt(nsubs(nd));
            yy(1,nd)=y(1,nd); ee(1,nd)=e(1,nd);
            labelcap{1} = 'FC';
            for c = 1:k
                y(c+1,nd) = FC_group_mean_at_cap{nd}(s1,s2,c);
                e(c+1,nd) = FC_group_std_at_cap{nd}(s1,s2,c)/sqrt(nsubs(nd));
                yy(c+1,nd) = FC_group_mean_regc{nd}(s1,s2,c);
                ee(c+1,nd) = FC_group_std_regc{nd}(s1,s2,c)/sqrt(nsubs(nd));
            end
            labelcap{c+1} = ['FC at C' num2str(c)];
        end
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 20+k;
        barwitherr(e,y) 
        xticks(1:k+2)
        xticklabels(labelcap)
        h_legend = legend(dslabel);
        set(h_legend,'Location','bestoutside')
        ylabel({'FC at CAPs (mean+/-SEM)'})
        box off
        print(['FC_at_cap_' s1name '-' s2name], '-dpng','-r600')
        clear x y e 
        close all
        
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 20+k;
        barwitherr(ee,yy) 
        xticks(1:k+2)
        xticklabels(labelcap)
        h_legend = legend(dslabel);
        set(h_legend,'Location','bestoutside')
        ylabel({'FC reg CAPs (mean+/-SEM)'})
        box off
        print(['FC_reg_cap_' s1name '-' s2name], '-dpng','-r600')
        clear x y e 
        close all
    end
end
    

end % function