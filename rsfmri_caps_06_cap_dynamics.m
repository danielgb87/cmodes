function caps_06_dFC_caps_compute_dynamics

% this function computes the CAP's timecourse and GS for a given k; their
% respective PSD; GS phase at the occurrence of a CAP; and performs
% comparisons between grousp. REEMEMBER THAT THE ORDER OF THE CAPS HERE IS
% THE SAME AS THE RESULTS FILE, AND AS THEY AREMAPPED. BUT THE ORDER OF THE
% PRELIMINARY "CAPS_PROPERTIES" AND "CAPS_QUALITY" ARE DIFFERENT AND
% ARRANGED ACCOORDING TO THE CAP OCCURRENCE RATE FROM THE CONCATENATED
% ANALYSIS. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 01-27-2021
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% path to dFC_CAPs_scripts and external toolboxes
addpath(genpath('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127'))

% 1.a load data, inputs, mask, and parcellation info.
main = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221';
main_analysis = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/analysis_vw_plus/cap_maps_lowres_k6';

addpath(genpath(main))
addpath(genpath(main_analysis))

cd(main_analysis)
k=6;
load(['cap_results_k' num2str(k)])
load('inputs')

% define and/or change datasets to include, names, and order of CAPs
% (depending on K)
r.ndatasets = length(inputs.ds_name); % number of datasets (session, groups, etc)
r.ds_name = inputs.ds_name;
r.rm_subjects = inputs.rm_subjects;
r.cap_ord=[3,6,2,5,1,4];
r.TR=inputs.analysis_TR;

% load the data to compute VE from, and rename it
cd(main)
load('data_vw_awk_lowres.mat')
data = data_vw;
clear data_vw;
if inputs.pproc_scrub==1
    clear data
    load('data_vw_scrub_awk_lowres.mat')
    data = data_vw_scrub;
    clear data_vw_scrub;

end
load('motion_info_0750_awk_lowres.mat')

ss=0;
for nd = 1:r.ndatasets
    r.nsubs(nd) = length(data{nd});
    if isempty(r.rm_subjects{nd}) == 0
        data{nd}(r.rm_subjects{nd}) = [];
        r.nsubs(nd)=r.nsubs(nd)-length(r.rm_subjects{nd});
    end
    for sub = 1:r.nsubs(nd)
        ss=ss+1;
        r.sub_list_ds{nd}(sub) = ss;
        for c = 1:k
            r.cap_occ{nd}(sub,c)= results_map.occ_prob_sub(ss,r.cap_ord(c));
            % put cap frames in the assigned order
            for f = 1:length(results_map.frame_ind_sub{ss})
                r.frame_ind_sub{nd}{sub}(f,1)=r.cap_ord(results_map.frame_ind_sub{ss}(f));
            end
        end
        r.btw_cap_sim{nd}(sub,:,:)=corr(results_map.map_vw_ds{nd}.cap_mean_map(r.cap_ord,:)');
    end
end
    
% finally choose filter band and cap timecourse threshold and name the
% folder to save results based on it (recommended)
r.f_low = 0.01;  % low bound frequency of filter for GS-phase
r.f_up = 0.03;   % upper bound frequency of filter for GS-phase
r.FS = 1./r.TR;  % sampling frequency for PSD
ffmax=0.15;      %maximum frequency in PSD plots
r.cfc_thr = 1;   % cap-to-frame correaltion thrshold above which to sample a CAP occurrences
r.xcorr_TT=30;   % max lag for cross-correlations
cd(main_analysis)
mkdir(['caps_dynamics_paper_k' num2str(k) '_thr1_01-03hz'])
cd(['caps_dynamics_paper_k' num2str(k) '_thr1_01-03hz'])

%% 2. compute cap timecourses, occurrence rates, GS, and their psd

%% 2.1 CAP cfc, btw CAP similarity (spatial and temporal) and GS corr
for nd = 1:r.ndatasets
    [r.cfc{nd}, r.btw_cap_corr_temp{nd}] = dFC_utils_caps_cfc(data{nd}, results_map.map_vw_ds{nd}.cap_mean_map(r.cap_ord,:));
end

%compute GS and CAP to GS temporal correlations
for nd = 1:r.ndatasets
    for sub = 1:r.nsubs(nd)
        r.gs{nd}{sub,1} = mean(data{nd}{sub},2);
        r.gs_norm{nd} = dFC_pproc_normalize_data(r.gs{nd});
        for c = 1:k
            r.cap2gs_corr{nd}(sub,c)=corr(spm_vec(r.gs_norm{nd}{sub}),spm_vec(r.cfc{nd}.cfc_norm{sub,c}));
        end
    end
    for c =1:k
        r.cap2gs_corr_mean{nd}(c)=mean(r.cap2gs_corr{nd}(:,c));   
        r.cap2gs_corr_sem{nd}(c)=std(r.cap2gs_corr{nd}(:,c))/sqrt(r.nsubs(nd));
    end
    writematrix(r.cap2gs_corr{nd}, ['cap2gs_corr_' r.ds_name{nd} '.xlsx'])
end

%concat between CAP similarity
dFC_utils_plot_imagesc_matrix_with_values(results_map.btw_cap_sim(r.cap_ord,r.cap_ord), 'jet')
print(['btw_cap_sim_concat'], '-dpng','-r600')
close all
writematrix(results_map.btw_cap_sim(r.cap_ord,r.cap_ord), ['btw_cap_sim_concat.xlsx'])

for nd = 1:r.ndatasets
    dFC_utils_plot_imagesc_matrix_with_values(r.btw_cap_corr_temp{nd}, 'jet')
    print(['btw_cap_corr_temp_' r.ds_name{nd}], '-dpng','-r600')
    close all
    writematrix(r.btw_cap_corr_temp{nd}, ['btw_cap_corr_temp_' r.ds_name{nd} '.xlsx'])

    dFC_utils_plot_imagesc_matrix_with_values(squeeze(mean(r.btw_cap_sim{nd},1)), 'jet')
    print(['btw_cap_sim_spat_' r.ds_name{nd}], '-dpng','-r600')
    close all
    writematrix(squeeze(mean(r.btw_cap_sim{nd},1)), ['btw_cap_corr_temp_' r.ds_name{nd} '.xlsx'])

    
    figure
    bar(1:k, r.cap2gs_corr_mean{nd});
    hold on
    errorbar(1:k, r.cap2gs_corr_mean{nd},r.cap2gs_corr_sem{nd},'.')
    xlabel('CAP','fontsize',14)
    ylabel('CAP to GS corr (mean+/- SEM)','fontsize',14)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 14);
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 14);
    xlim([0 k+1])
    print(['CAP2GS_corr_' r.ds_name{nd}], '-dpng','-r600')
    close all
    writematrix(squeeze(r.cap2gs_corr_mean{nd}), ['CAP2GS_corr_' r.ds_name{nd} '.xlsx'])

end


% compute auto and cross-correlations
for nd=1:r.ndatasets
    for c1 = 1:k
        for c2 = 1:k
            for sub = 1:r.nsubs(nd)
                %xcorr
                r.cfc_xcorr{nd}{c1,c2}(sub,:) = xcorr(spm_vec(r.cfc{nd}.cfc_norm{sub,c1}),spm_vec(r.cfc{nd}.cfc_norm{sub,c2}),r.xcorr_TT, 'coeff');
            end
            r.cfc_xcorr_mean{nd}{c1,c2}(:)=mean(r.cfc_xcorr{nd}{c1,c2},1);
            r.cfc_xcorr_sem{nd}{c1,c2}(:)=std(r.cfc_xcorr{nd}{c1,c2},1)/sqrt(r.nsubs(nd));
           % writematrix(r.cfc_xcorr{nd}{c1,c2}, ['xcorr_' r.ds_name{nd} '_cap' num2str(c1) 'to' num2str(c2) '.xlsx'])
            
%             figure
%             fig = gcf;
%             fig.Units = 'centimeters';
%             fig.Position(3) = 13;
%             fig.Position(4) = 6;
%             fig.Position(1) = 20;
%             fig.Position(2) = 20;
%             if c1 == c2
%                 shadedErrorBar(-r.xcorr_TT*r.TR(nd):r.TR(nd):r.xcorr_TT*r.TR(nd),r.cfc_xcorr_mean{nd}{c1,c2},r.cfc_xcorr_sem{nd}{c1,c2},'r')
%             else
%                 shadedErrorBar(-r.xcorr_TT*r.TR(nd):r.TR(nd):r.xcorr_TT*r.TR(nd),r.cfc_xcorr_mean{nd}{c1,c2},r.cfc_xcorr_sem{nd}{c1,c2},'b')
%             end
%             xlim([-r.xcorr_TT*r.TR(nd) ,r.xcorr_TT*r.TR(nd)])
%             a=get(gca,'XTickLabel');
%             set(gca,'XTickLabel',a,'fontsize',10)
%             ylim([-1,1])
%             xlabel('lag (s)', 'FontSize', 10)
%             ylabel(['cross-correlation (mean +/-SEM)'], 'FontSize', 10)
%             box on, grid on
%             fig.PaperPositionMode = 'auto';
%             print(['xcorr_' r.ds_name{nd} '_cap' num2str(c1) 'to' num2str(c2)], '-dpng','-r600')
%             close all

        end
        
    end
   
end



%% 2.2 cap occurrence rate
%CAP occurrence rate
for nd = 1:r.ndatasets
    for c = 1:k
        r.cap_occ_mean{nd}(c) = mean(r.cap_occ{nd}(:,c));
        r.cap_occ_sem{nd}(c) = std(r.cap_occ{nd}(:,c))/sqrt(r.nsubs(nd));
        label{c}=['c' num2str(c)];
    end
    writematrix(squeeze(r.cap_occ{nd}), ['cap_occ_' r.ds_name{nd} '.xlsx'])
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 13;
    fig.Position(4) = 10; fig.Position(1) = 20; fig.Position(2) = 20;
    bar(r.cap_occ_mean{nd})
    hold on
    errorbar(1:k,r.cap_occ_mean{nd},r.cap_occ_sem{nd},'.')
    xlim([0 k+1])
    set(gca,'XTickLabel',label)
    ylabel({'occurrence rate (mean+/- SEM)'})
    box off
    fig.PaperPositionMode = 'auto';
    print(['cap_occ_' r.ds_name{nd}], '-dpng','-r600')
    close all
end

%% 2.3 r.motion2CAP distributions
% optional, if motion was not scrubbed, compute the proportion of motion
% frames in a CAP.
if inputs.pproc_scrub == 0
    for nd = 1:r.ndatasets
        
        [r.motion2CAP{nd}] = dFC_utils_motion2cap(motion_info{nd}, results_map.frame_ind_sub(r.sub_list_ds{nd}), r.rm_subjects{nd});    
        if isempty(r.rm_subjects{nd}) == 0
            motion_info.FD_flag(r.rm_subjects{nd}) = [];
            motion_info.FD_flag_proportion(r.rm_subjects{nd}) = [];
            motion_info.FD(r.rm_subjects{nd}) = [];
        end
        for sub = 1:r.nsubs(nd)
            r.motion2CAP{nd}.count{sub}= r.motion2CAP{nd}.count{sub}(:,r.cap_ord);
            r.motion2CAP{nd}.proportion_allTS= r.motion2CAP{nd}.proportion_allTS(:,r.cap_ord);
            r.motion2CAP{nd}.cap_proportion= r.motion2CAP{nd}.cap_proportion(:,r.cap_ord);
            r.motion2CAP{nd}.count_all= r.motion2CAP{nd}.count_all(:,r.cap_ord);
            r.motion2CAP{nd}.FD_flag_at_cap{sub}= r.motion2CAP{nd}.FD_flag_at_cap{sub}(:,r.cap_ord);
            r.motion2CAP{nd}.FD_flag_at_cap_mean= r.motion2CAP{nd}.FD_flag_at_cap_mean(:,r.cap_ord);
            r.motion2CAP{nd}.FD_at_cap{sub}= r.motion2CAP{nd}.FD_at_cap{sub}(:,r.cap_ord);
            r.motion2CAP{nd}.FD_at_cap_mean= r.motion2CAP{nd}.FD_at_cap_mean(:,r.cap_ord);
        end
        writematrix(squeeze(r.motion2CAP{nd}.FD_at_cap_mean), ['FD_at_cap_mean_' r.ds_name{nd} '.xlsx'])
        writematrix(squeeze(r.motion2CAP{nd}.FD_flag_at_cap_mean), ['FD_flag_at_cap_mean_' r.ds_name{nd} '.xlsx'])
        writematrix(squeeze(r.motion2CAP{nd}.cap_proportion), ['motion2CAP_norm_' r.ds_name{nd} '.xlsx'])
        writematrix(squeeze(r.motion2CAP{nd}.proportion_allTS), ['motion2CAP_' r.ds_name{nd} '.xlsx'])

        
        [r.motion2CAP{nd}.motion2CAP_anova_p,r.motion2CAP{nd}.motion2CAP_anova_table,s1] = anova1(r.motion2CAP{nd}.proportion_allTS);
        [r.motion2CAP{nd}.motion2CAP_proportion_anova_p,r.motion2CAP{nd}.motion2CAP_proportion_anova_table,s2] = anova1(r.motion2CAP{nd}.cap_proportion);
        [r.motion2CAP{nd}.motion2CAP_FD_anova_p,r.motion2CAP{nd}.motion2CAP_FD_anova_table,s3] = anova1(r.motion2CAP{nd}.FD_at_cap_mean);
        close all
        [r.motion2CAP{nd}.motion2CAP_anova_multComp] = multcompare(s1,'alpha',0.05,'ctype','hsd');
        saveas(gcf,['CAP2motion_anova_multComp' r.ds_name{nd} '.png'])
        [r.motion2CAP{nd}.motion2CAP_proportion_anova_multComp] = multcompare(s2,'alpha',0.05,'ctype','hsd');
        saveas(gcf,['CAP2motion_proportion_anova_multComp' r.ds_name{nd} '.png'])
        [r.motion2CAP{nd}.motion2CAP_FD_anova_multComp] = multcompare(s3,'alpha',0.05,'ctype','hsd');
        saveas(gcf,['CAP2motion_FD_anova_multComp' r.ds_name{nd} '.png'])
        close all

    end
end


%% 2.4 CAP and GS PSD
% compute psd with fourier and multitaper methods
for nd = 1:r.ndatasets
    r.cfc_psd{nd} = dFC_utils_caps_psd(r.cfc{nd}.cfc_norm, r.TR(nd));
    r.gs_psd{nd} = dFC_utils_caps_psd(r.gs_norm{nd}, r.TR(nd));
    r.NW = 4;
    r.cfc_psdMT{nd} = dFC_utils_caps_psdMT(r.cfc{nd}.cfc_norm,r.TR(nd),r.NW);
    r.gs_psdMT{nd} = dFC_utils_caps_psdMT(r.gs_norm{nd}, r.TR(nd),r.NW);
    
    for c1 = 1:k
        figure,
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
        fig.Position(4) = 6; fig.Position(1) = 20; fig.Position(2) = 20;
        shadedErrorBar(r.cfc_psdMT{nd}.min_freq_vec, r.cfc_psdMT{nd}.psd_ds_mean{c1},r.cfc_psdMT{nd}.psd_ds_std{c1}/sqrt(r.nsubs(nd)),'b')
        [M,ii] =max(r.cfc_psdMT{nd}.psd_ds_mean{c1}); fp = r.cfc_psdMT{nd}.min_freq_vec(ii);
        ymax = max(max(spm_vec(r.cfc_psdMT{nd}.psd_ds_mean)))+1.2*max(max(spm_vec(r.cfc_psdMT{nd}.psd_ds_std)))/sqrt(r.nsubs(nd));
        ylim([0,round(ymax,-1)]), xlim([0 ffmax])
        text(double(fp*1.1), 1.1*double(M), [num2str(round(fp,3)) ' Hz'],'FontSize', 10);
        % add vertical lines for the 0.01-0.03 Hz band
        hold on
        plot(0.01*ones(round(ymax,-1),1),1/round(ymax,-1):round(ymax,-1)/round(ymax,-1):round(ymax,-1),'--r')
        hold on
        plot(0.03*ones(round(ymax,-1),1),1/round(ymax,-1):round(ymax,-1)/round(ymax,-1):round(ymax,-1),'--r')
        box off
        legend({['CAP ' num2str(c1)]},'FontSize', 10)
        xlabel('Frequency (Hz)', 'FontSize', 10)
        ylabel(['PSD (mean +/-SEM)'], 'FontSize', 10)
        box off
        fig.PaperPositionMode = 'auto';
        print(['psdMT_' r.ds_name{nd} '_cap' num2str(c1)], '-dpng','-r600')
        close all
    end
    % now GS PSD
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
    fig.Position(4) = 6; fig.Position(1) = 20; fig.Position(2) = 20;
    shadedErrorBar(r.gs_psdMT{nd}.min_freq_vec, r.gs_psdMT{nd}.psd_ds_mean{1},r.gs_psdMT{nd}.psd_ds_std{1}/sqrt(r.nsubs(nd)),'b')
    [M,ii] =max(r.gs_psdMT{nd}.psd_ds_mean{1}); fp = r.gs_psdMT{nd}.min_freq_vec(ii);
    ymax = max(max(spm_vec(r.gs_psdMT{nd}.psd_ds_mean{1})))+1.2*max(max(spm_vec(r.gs_psdMT{nd}.psd_ds_std{1})))/sqrt(r.nsubs(nd));
    ylim([0,round(ymax,-1)]), xlim([0 ffmax])
    text(double(fp*1.1), 1.1*double(M), [num2str(round(fp,3)) ' Hz'],'FontSize', 10);
    % add vertical lines for the 0.01-0.03 Hz band
    hold on
    plot(0.01*ones(round(ymax,-1),1),1/round(ymax,-1):round(ymax,-1)/round(ymax,-1):round(ymax,-1),'--r')
    hold on
    plot(0.03*ones(round(ymax,-1),1),1/round(ymax,-1):round(ymax,-1)/round(ymax,-1):round(ymax,-1),'--r')
    box off
    legend({['GS']},'FontSize', 10)
    xlabel('Frequency (Hz)', 'FontSize', 10)
    ylabel(['PSD (mean +/-SEM)'], 'FontSize', 10)
    box off
    fig.PaperPositionMode = 'auto';
    print(['psdMT_' r.ds_name{nd} '_gs'], '-dpng','-r600')
    close all
    
    % save as xlsx format
    for sub = 1:r.nsubs(nd)
        psd_gs_mat{nd}(:,1)=r.gs_psdMT{nd}.min_freq_vec;
        for c = 1:k
            psd_mat{nd,c}(:,1)=r.cfc_psdMT{nd}.min_freq_vec;
            psd_mat{nd,c}(:,1+sub)=r.cfc_psdMT{nd}.psd_ds{sub,c};
        end
        psd_gs_mat{nd}(:,1+sub)=r.gs_psdMT{nd}.psd_ds{sub};
    end
    for c = 1:k
        writematrix(squeeze(psd_mat{nd,c}), ['psdMT_' r.ds_name{nd} '_cap' num2str(c) '.xlsx'])
    end
    writematrix(squeeze(psd_gs_mat{nd}), ['psdMT_GS_' r.ds_name{nd} '.xlsx'])

    
end


%% 3. filter GS and PSD and compute the GS phase distributions at the occurrence of each CAP
% filter GS and CFC and compute the GS phases at cap...also do it after
% concatenating all datasets
ss=0;
for nd = 1:r.ndatasets

    [r.gs_filt{nd}] = dFC_utils_caps_NBfilter(r.gs_norm{nd},r.f_low,r.f_up,r.FS(nd));
    [r.cfc_filt{nd}] = dFC_utils_caps_NBfilter(r.cfc{nd}.cfc_norm,r.f_low,r.f_up,r.FS(nd));
    [r.gs_phase{nd}] = dFC_utils_caps_comp_phase(r.gs_filt{nd});
    [r.gs_phase_at_cap{nd}] = dFC_utils_caps_gs_phase_at_cap(r.gs_phase{nd}.phase, r.cfc_filt{nd}, r.frame_ind_sub{nd}, r.cfc_thr);
    [r.gs_phase_at_cap_wtn_cycles{nd}] = dFC_utils_gs_phase_at_cap_in_cycle(r.gs_phase{nd}.phase, r.cfc_filt{nd},r.frame_ind_sub{nd}, r.cfc_thr);
    
    for sub = 1:r.nsubs(nd)
        ss = ss+1;
        r.gs_concat{ss,1} = r.gs{nd}{sub};
        r.gs_norm_concat{ss,1} = r.gs_norm{nd}{sub};
        r.gs_filt_concat{ss,1} = r.gs_filt{nd}{sub};
        r.gs_phase_concat{ss,1} = r.gs_phase{nd}.phase{sub};
        r.frame_ind_sub_concat{ss}=r.frame_ind_sub{nd}{sub};
        for c = 1:k
            r.cfc_concat{ss,c} = r.cfc{nd}.cfc{sub,c};
            r.cfc_norm_concat{ss,c} = r.cfc{nd}.cfc_norm{sub,c};
            r.cfc_filt_concat{ss,c} = r.cfc_filt{nd}{sub,c};
        end
    end
    
    for c = 1:k
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 5;
        fig.Position(4) = 5; fig.Position(1) = 20; fig.Position(2) = 20;
        polarhistogram(spm_vec(r.gs_phase_at_cap_wtn_cycles{nd}.gs_phase_at_CAP_thr{c}),20,'normalization','probability')
        pax = gca; pax.ThetaTick = [0 90 180 270] ; pax.FontSize = 12; pax.RTick = [0.05 0.1 0.15];
        pax.ThetaAxisUnits = 'radians'; pax.ThetaColor = 'b'; pax.RColor = 'k';
        pax.LineWidth = 2;
        title(['GS phase at CAP ' num2str(c)], 'FontSize',12)
        set(gcf,'color','none')
        set(gca,'color','none')
        fig.PaperPositionMode = 'auto';
        print(['gs_phase_' r.ds_name{nd} '_cap' num2str(c)], '-dpng','-r600')
        close all
        label{c} = ['C' num2str(c)];
    end
end

[r.gs_phase_at_cap_concat] = dFC_utils_caps_gs_phase_at_cap(r.gs_phase_concat, r.cfc_filt_concat, r.frame_ind_sub_concat, r.cfc_thr);
[r.gs_phase_at_cap_wtn_cycles_concat] = dFC_utils_gs_phase_at_cap_in_cycle(r.gs_phase_concat, r.cfc_filt_concat,r.frame_ind_sub_concat, r.cfc_thr);
for c = 1:k
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 5;
    fig.Position(4) = 5; fig.Position(1) = 20; fig.Position(2) = 20;
    polarhistogram(spm_vec(r.gs_phase_at_cap_wtn_cycles_concat.gs_phase_at_CAP_thr{c}),20,'normalization','probability')
    pax = gca; pax.ThetaTick = [0 90 180 270] ; pax.FontSize = 12; pax.RTick = [0.05 0.1 0.15];
    pax.ThetaAxisUnits = 'radians'; pax.ThetaColor = 'b'; pax.RColor = 'k';
    pax.LineWidth = 2;
    title(['GS phase at CAP ' num2str(c)], 'FontSize',12)
    set(gcf,'color','none')
    set(gca,'color','none')
    fig.PaperPositionMode = 'auto';
    print(['gs_phase_CONCAT_cap' num2str(c)], '-dpng','-r600')
    close all
end


%% 4. compute fALFF
% !!!!!!!!!!! BY DEFAULT IT IS COMPUTE IN THE 0.01-0.03 ti 0.01-0.1 hz
% ratio. IF YOUR POWER SPECTRA DIFFERS, THEN EDIT THESE VALUES BELOW.
for nd = 1:r.ndatasets
    for c = 1:k 
        label{c} = ['C' num2str(c)];

        %cap and GS ALFF.
        for sub = 1:r.nsubs(nd)
            r.cap_falff{nd}(c,sub) = bandpower(r.cfc_psdMT{nd}.psd{sub,c},r.cfc_psdMT{nd}.freqs{sub},[0.01 0.03],'psd')/bandpower(r.cfc_psdMT{nd}.psd{sub,c},r.cfc_psdMT{nd}.freqs{sub},[0.01 0.1],'psd');       
            r.gs_falff{nd}(sub) = bandpower(r.gs_psdMT{nd}.psd{sub},r.gs_psdMT{nd}.freqs{sub},[0.01 0.03],'psd')/bandpower(r.gs_psdMT{nd}.psd{sub},r.gs_psdMT{nd}.freqs{sub},[0.01 0.1],'psd');
        end
        r.cap_falff_mean{nd}(c)=mean(r.cap_falff{nd}(c,:));
        r.cap_falff_std{nd}(c)=std(r.cap_falff{nd}(c,:));
        r.gs_falff_mean(nd)=mean(r.gs_falff{nd}(:));
        r.gs_falff_std(nd)=std(r.gs_falff{nd}(:));
    end
    writematrix(squeeze(r.cap_falff{nd}), ['falff_caps_' r.ds_name{nd} '.xlsx'])
    writematrix(squeeze(r.gs_falff{nd}), ['falff_gs_' r.ds_name{nd} '.xlsx'])

    r.falff_mean{nd} = [r.cap_falff_mean{nd} r.gs_falff_mean(nd)];
    r.falff_std{nd} = [r.cap_falff_std{nd} r.gs_falff_std(nd)];
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 13;
    fig.Position(4) = 10; fig.Position(1) = 20; fig.Position(2) = 20;
    bar(r.falff_mean{nd})
    hold on
    errorbar(1:k+1,r.falff_mean{nd},r.falff_std{nd}/sqrt(r.nsubs(nd)),'.')
    xlim([0 k+2])
    set(gca,'XTickLabel',label)
    ylabel({'CAP and GS fALFF (mean+/- SEM)'})
    box off
    fig.PaperPositionMode = 'auto';
    print(['fALFF_' r.ds_name{nd}], '-dpng','-r600')
    close all
end



%% 5. ORGANIZE INFORMATION IN A STRUCTURE
caps_dyn = r;
save('caps_dyn','caps_dyn')



end % function