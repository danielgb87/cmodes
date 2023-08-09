function dFC_caps_compute_dynamics

% this function computes the CAP's timecourse and GS for a given k; their
% respective PSD; GS phase at the occurrence of a CAP; and performs
% comparisons between grousp. REEMEMBER THAT THE ORDER OF THE CAPS HERE IS
% THE SAME AS THE RESULTS FILE, AND AS THEY AREMAPPED. BUT THE ORDER OF THE
% PRELIMINARY "CAPS_PROPERTIES" AND "CAPS_QUALITY" ARE DIFFERENT AND
% ARRANGED ACCOORDING TO THE CAP OCCURRENCE RATE FROM THE CONCATENATED
% ANALYSIS. 

clc,clear
%% EDIT HERE
%% 1.a load caps results, inputs, mask, and parcellation info.
k=12;
main = ['/media/DATA/dgutierrez/human_data/caps_hnu_10session_GM_concat_3mm/analysis_parcel/cap_maps_3mm_k' num2str(k)];
cd(main)
ndatasets = 10; % define the amount of datasets to compare
ncaps = 2:20;    % the clusterings you want to compare... make sure all datasets and runs have them.
main_ds = '/media/DATA/dgutierrez/human_data/caps_hnu_10session_GM_concat_3mm';
dataset_name = 'concat_parcel';
nruns = 5;
%Load the results from each datasets for the best run
cd(main)
load(['cap_results_k' num2str(k)])
load('inputs')
rm_subjects = inputs.rm_subjects;

% load the data to compute VE from, and rename it
cd(main_ds)
load('data_vw_hnu_mni152_3mm_GM.mat')
load('motion_info_hnu_mni152_3mm_GM.mat')

TR=inputs.analysis_TR;
data = data_vw;
clear data_vw;
ss=0;
for nd = 1:ndatasets
    nsubs(nd) = length(data{nd});
    if isempty(rm_subjects{nd}) == 0
        data{nd}(rm_subjects{nd}) = [];
        nsubs(nd)=nsubs(nd)-length(rm_subjects{nd});
    end
    for sub = 1:nsubs(nd)
        ss=ss+1;
        sub_list_ds{nd}(sub) = ss;
    end
end

% finally choose filter band and cap timecourse threshold 
f_low = 0.06;
f_up = 0.1;
FS = 1/TR;
cfc_thr = 1;
cd(main)
mkdir(['caps_dynamics_k' num2str(k) '_thr1_06-1hz'])
cd(['caps_dynamics_k' num2str(k) '_thr1_06-1hz'])

%% 2. compute cap timecourses, GS, and their psd

for nd = 1:ndatasets
    [cfc{nd}, btw_cap_corr{nd}] = dFC_utils_caps_cfc(data{nd}, results_map.map_vw_ds{nd}.cap_mean_map);
end

% optional, if motion was not scrubbed, compute the proportion of motion
% frames in a CAP.
if inputs.pproc_scrub == 0
    for nd = 1:ndatasets
        [motion2CAP{nd}] = dFC_utils_motion2cap(motion_info{nd}, results_map.map_vw.frame_ind_sub(sub_list_ds{nd}), rm_subjects{nd});
    end
end

%compute GS
for nd = 1:ndatasets
    for sub = 1:nsubs(nd)
        gs{nd}{sub,1} = mean(data{nd}{sub},2);
        gs_norm{nd} = dFC_pproc_normalize_data(gs{nd});
    end
end

% compute psd with fourier and multitaper methods
for nd = 1:ndatasets
    cfc_psd{nd} = dFC_utils_caps_psd(cfc{nd}.cfc_norm, TR);
    gs_psd{nd} = dFC_utils_caps_psd(gs_norm{nd}, TR);
    NW = 4;
    cfc_psdMT{nd} = dFC_utils_caps_psdMT(cfc{nd}.cfc_norm,TR,NW);
    gs_psdMT{nd} = dFC_utils_caps_psdMT(gs_norm{nd}, TR,NW);
end

% plot PSD
% CAPs
for nd = 1:ndatasets
    for c1 = 1:k
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
        fig.Position(4) = 6; fig.Position(1) = 20; fig.Position(2) = 20;   
        shadedErrorBar(cfc_psdMT{nd}.min_freq_vec, cfc_psdMT{nd}.psd_ds_mean{c1},cfc_psdMT{nd}.psd_ds_std{c1}/sqrt(nsubs(nd)),'b')
        ([M,ii] =max(cfc_psdMT{nd}.psd_ds_mean{c1}); fp = cfc_psdMT{nd}.min_freq_vec(ii);
        ymax = max(max(spm_vec(cfc_psdMT{nd}.psd_ds_mean)))+1.2*max(max(spm_vec(cfc_psdMT{nd}.psd_ds_std)))/sqrt(nsubs(nd));
        ylim([0,round(ymax,-1)]), xlim([0 0.15])
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
        print(['psdMT_ds' num2str(nd) '_cap' num2str(c1)], '-dpng','-r600')
        close all
    end    
    % now GS PSD
    figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
        fig.Position(4) = 6; fig.Position(1) = 20; fig.Position(2) = 20;   
        shadedErrorBar(gs_psdMT{nd}.min_freq_vec, gs_psdMT{nd}.psd_ds_mean{1},gs_psdMT{nd}.psd_ds_std{1}/sqrt(nsubs(nd)),'b')
        [M,ii] =max(gs_psdMT{nd}.psd_ds_mean{1}); fp = gs_psdMT{nd}.min_freq_vec(ii);
        ymax = max(max(spm_vec(gs_psdMT{nd}.psd_ds_mean{1})))+1.2*max(max(spm_vec(gs_psdMT{nd}.psd_ds_std{1})))/sqrt(nsubs(nd));
        ylim([0,round(ymax,-1)]), xlim([0 0.15])
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
        print(['psdMT_ds' num2str(nd) '_gs'], '-dpng','-r600')
        close all
end

% plot between CAP corr
for nd = 1:ndatasets
    dFC_utils_plot_imagesc_matrix_with_values(btw_cap_corr{nd}, 'jet')
    print(['btw_cap_corr_ds' num2str(nd)], '-dpng','-r600')
    close all
end

%plot motion2CAP
for nd=1:ndatasets
    if isempty(rm_subjects) == 0
    motion_info{nd}.FD_flag_proportion(rm_subjects{nd}) = [];
    end
    s = 0;
    flag_vols = motion_info{nd}.FD_flag_proportion;
    [fl subj] = sort(flag_vols, 'descend');
    fl = fl(1:nsubs(nd));
    for sub = 1:length(fl)
        motion_cap{nd}(sub,:) = motion2CAP{nd}.count_all(subj(sub),:);
    end
    [motion2CAP{nd}.m2cap_anova_table,a,s] = anova1(motion_cap{nd});
    close all
    [motion2CAP{nd}.m2cap_anova_multComp] = multcompare(s,'alpha',0.05,'ctype','hsd');
    saveas(gcf,['CAP2motion_multComp' num2str(nd) '.fig'])
    close all
    figure
    bar(1:k, mean(motion_cap{nd},1));
    hold on
    errorbar(1:k, mean(motion_cap{nd},1),std(motion_cap{nd},1)/sqrt(nsubs(nd)),'.')
    xlabel('CAP','fontsize',14)
    ylabel('Motion volumes per CAP (+/- SEM)','fontsize',14)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 14);
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 14);
    xlim([0 k+1])
    print(['CAP2motion_' num2str(nd)], '-dpng','-r600')
    close all
end


%% 3. filter GS and PSD and compute the GS phase distributions at the occurrence of each CAP
% filter GS and CFC

for nd = 1:ndatasets
    [gs_filt{nd}] = dFC_utils_caps_NBfilter(gs_norm{nd},f_low,f_up,FS);
    [cfc_filt{nd}] = dFC_utils_caps_NBfilter(cfc{nd}.cfc_norm,f_low,f_up,FS);
    [gs_phase{nd}] = dFC_utils_caps_comp_phase(gs_filt{nd});
    [gs_phase_at_cap{nd}] = dFC_utils_caps_gs_phase_at_cap(gs_phase{nd}.phase, cfc_filt{nd}, results_map.frame_ind_sub(sub_list_ds{nd}), cfc_thr);
    [gs_phase_at_cap_wtn_cycles{nd}] = dFC_utils_gs_phase_at_cap_in_cycle(gs_phase{nd}.phase, cfc_filt{nd},results_map.frame_ind_sub(sub_list_ds{nd}), cfc_thr);
end

%% 4. compute fALFF
% !!!!!!!!!!! BY DEFAULT IT IS COMPUTE IN THE 0.01-0.03 ti 0.01-0.1 hz
% ratio. IF YOUR POWER SPECTRA DIFFERS, THEN EDIT THESE VALUES BELOW.
for nd = 1:ndatasets
    for c = 1:k 
        %cap and GS ALFF.
        for sub = 1:nsubs(nd)
            cap_falff{nd}(c,sub) = bandpower(cfc_psdMT{nd}.psd{sub,c},cfc_psdMT{nd}.freqs{sub},[0.01 0.03],'psd')/bandpower(cfc_psdMT{nd}.psd{sub,c},cfc_psdMT{nd}.freqs{sub},[0.01 0.1],'psd');       
            gs_falff{nd}(sub) = bandpower(gs_psdMT{nd}.psd{sub},gs_psdMT{nd}.freqs{sub},[0.01 0.03],'psd')/bandpower(gs_psdMT{nd}.psd{sub},gs_psdMT{nd}.freqs{sub},[0.01 0.1],'psd');
        end
        cap_falff_mean{nd}(c)=mean(cap_falff{nd}(c,:));
        cap_falff_std{nd}(c)=std(cap_falff{nd}(c,:));
        gs_falff_mean(nd)=mean(gs_falff{nd}(:));
        gs_falff_std(nd)=std(gs_falff{nd}(:));
    end
    mean_falff{nd} = [cap_falff_mean{nd} gs_falff_mean(nd)];
    std_falff{nd} = [cap_falff_std{nd} gs_falff_std(nd)];
end

%plot GS phase distributions and fALFF
for nd = 1:ndatasets
    for c = 1:k
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 5;
        fig.Position(4) = 5; fig.Position(1) = 20; fig.Position(2) = 20;
        polarhistogram(spm_vec(gs_phase_at_cap_wtn_cycles{nd}.gs_phase_at_CAP_thr{c}),20,'normalization','probability')
        pax = gca; pax.ThetaTick = [0 90 180 270] ; pax.FontSize = 12; pax.RTick = [0.05 0.1 0.15];
        pax.ThetaAxisUnits = 'radians'; pax.ThetaColor = 'b'; pax.RColor = 'k';
        pax.LineWidth = 2;
        title(['GS phase at CAP ' num2str(c)], 'FontSize',12)
        set(gcf,'color','none')
        set(gca,'color','none')
        fig.PaperPositionMode = 'auto';
        print(['gs_phase_at_ds' num2str(nd) '_cap' num2str(c)], '-dpng','-r600')
        close all
        label{c} = ['C' num2str(c)];
    end
    label{c+1} = 'GS';
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 13;
    fig.Position(4) = 10; fig.Position(1) = 20; fig.Position(2) = 20;
    bar(mean_falff{nd})
    hold on
    errorbar(1:k+1,mean_falff{nd},std_falff{nd}/sqrt(nsubs(nd)),'.')
    xlim([0 k+2])
    set(gca,'XTickLabel',label)
    ylabel({'CAP and GS fALFF (mean+/- SEM)'})
    box off
    fig.PaperPositionMode = 'auto';
    print(['fALFF_ds' num2str(nd)], '-dpng','-r600')
    close all
end

%% 5. ORGANIZE INFORMATION IN A STRUCTURE
caps_dyn.cfc.cfc = cfc;
caps_dyn.cfc.cfc_filt = cfc_filt;
caps_dyn.cfc.cap_falff = cap_falff;
caps_dyn.cfc.psd = cfc_psd;
caps_dyn.cfc.psdMT = cfc_psdMT;
caps_dyn.cfc.cfc_threshold = cfc_thr;
caps_dyn.TR = TR;
caps_dyn.filter.flow = f_low;
caps_dyn.filter.fup = f_up;
caps_dyn.nsubs = nsubs;
caps_dyn.gs.gs = gs;
caps_dyn.gs.gs_norm = gs_norm;
caps_dyn.gs.gs_filt = gs_filt;
caps_dyn.gs.gs_falff = gs_falff;
caps_dyn.gs.gs_phase = gs_phase;
caps_dyn.gs.gs_psd = gs_psd;
caps_dyn.gs.gs_psdMT = gs_psdMT;
caps_dyn.gs.gs_phase_at_cap = gs_phase_at_cap;
caps_dyn.gs.gs_phase_at_cap_wtn_cycle = gs_phase_at_cap_wtn_cycles;
caps_dyn.btw_cap_corr = btw_cap_corr;
caps_dyn.motion2CAP = motion2CAP;
save('caps_dyn','caps_dyn')


end % function