function dFC_SBA_hrf_fALFF_psd


% This function takes a group of seeds, extract the data from each dataset, and computes the PSD, fALFF, and HRF using a double-gamma canonical function set.
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
TR=inputs.analysis_TR;
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
    for sub=1:nsubs(nd)
            nobs{nd}(sub) = length(seed_data{nd,1}.TS{sub});
    end
end

nseeds = length(seed_list);
cd(main)
%% 2. compute PSD and fALFF
% finally, compute fALFF for the selected seeds and plot the power spectra
NW = 4;
for s = 1:nseeds
    for nd=1:ndatasets
        NFFT = max(nobs{nd}); 
        [psd_MT{nd,s}] = dFC_utils_caps_psdMT(seed_data{nd,s}.TS, TR, NW);
        for sub = 1:nsubs(nd)
            roi_pband{nd}(s,sub) = bandpower(psd_MT{nd,s}.psd{sub},psd_MT{nd,s}.freqs{sub},[0.01 0.03],'psd');
            roi_ptot{nd}(s,sub) = bandpower(psd_MT{nd,s}.psd{sub},psd_MT{nd,s}.freqs{sub},[0.01 0.1],'psd');
            roi_alff{nd}(s,sub) = roi_pband{nd}(s,sub)/roi_ptot{nd}(s,sub);
        end        
        %stats on falff
        roi_alff_mean(nd,s) = mean(roi_alff{nd}(s,:));
        roi_alff_std(nd,s) = std(roi_alff{nd}(s,:));
        
    end
    for n1 = 1:ndatasets  
        for n2 = 1:ndatasets
            [roi_alff_h(n1,n2,s), roi_alff_p(n1,n2,s),~ , stt] = ttest2(spm_vec(roi_alff{n1}(s,:)), spm_vec(roi_alff{n2}(s,:)),'tail','both');
            roi_alff_T(n1,n2,s) = stt.tstat; clear stt
        end
    end
end
% plot PSD and fALFF        
for s = 1:nseeds
    for nd = 1:ndatasets
        figure
        fig = gcf;
        fig.Units = 'centimeters';
        fig.Position(3) = 8;
        fig.Position(4) = 6;
        fig.Position(1) = 5;
        fig.Position(2) = 5;   
        % start making figure     
        shadedErrorBar(psd_MT{nd,s}.min_freq_vec, psd_MT{nd,s}.psd_ds_mean{1},psd_MT{nd,s}.psd_ds_std{1}/sqrt(nsubs(nd)))        
        %call max psd and its freq.
        [M, ii] =max(psd_MT{nd,s}.psd_ds_mean{1});
        fp  = psd_MT{nd,s}.min_freq_vec(ii);
        %put limits
        xlim([0 0.11])
        xlabel('Frequency (Hz)', 'FontSize', 10)
        ylabel(['PSD (mean +/-SEM)'], 'FontSize', 10)
        box off
        fig.PaperPositionMode = 'auto';
        print(['PSD_ds' num2str(nd) '_' seed_list(s).name(1:end-7)], '-dpng','-r600')
        close all
    end
end


%% 3. compute hrfs

min_onset_search = 2; % minimum delay allowed between event and HRF onset (seconds)
max_onset_search = 8; % maximum delay allowed between event and HRF onset (seconds)
temporal_mask = [];
para_canon.estimation = 'canon';
para_canon.TR = TR;
para_canon.T  = 5; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid
para_canon.T0=fix(para_canon.T/2); % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)
para_canon.dt  = para_canon.TR/para_canon.T; % fine scale time resolution.
para_canon.TD_DD = 2; % time and dispersion derivative
para_canon.AR_lag = 1; % AR(1) noise autocorrelation.
para_canon.thr = 1; % (mean+) para.thr*standard deviation threshold to detect event.
para_canon.len = 18; % length of HRF, in seconds
para_canon.lag  = fix(min_onset_search/para_canon.dt):fix(max_onset_search/para_canon.dt);

for r = 1:nseeds
    for nd = 1:ndatasets
        for sub = 1:nsubs(nd)
            [beta_hrf{nd}{sub,r}, bf{nd}{sub,r}, event_bold{nd}{sub,r}] = wgr_rshrf_estimation_canonhrf2dd_par2(spm_vec(seed_data{nd,r}.TS{sub}),para_canon,temporal_mask);
            hrfa{nd,r}(sub,:) = bf{nd}{sub,r}*beta_hrf{nd}{sub,r}(1:size(bf{nd}{sub,r},2),:); %HRF canonical with 2 derivatives
            [params_hrf_canon{nd,r}(sub,:)] = wgr_get_parameters(hrfa{nd,r}(sub,:),para_canon.TR/para_canon.T);% estimate HRF parameter
        end
    end
end

%obtain mean HRFs and compute stats on parameters
for r = 1:nseeds
    for nd =1:ndatasets
        % time-lock average th HRFs (mean and SEM)
        for tt = 1:length(hrfa{nd,r}(1,:))
            hrf_mean{nd,r}(tt) = mean(hrfa{nd,r}(:,tt));
            hrf_std{nd,r}(tt) = std(hrfa{nd,r}(:,tt));
        end
    end
end


%plot and save the mean hrfs
for r = 1:nseeds 
    for nd = 1:ndatasets
    figure
        hold on
        shadedErrorBar(0:TR/para_canon.T:para_canon.len,hrf_mean{nd,r},hrf_std{nd,r}/sqrt(nsubs(nd)),'k',0.25)
    
    %set(gcf, 'InvertHardCopy','off')
    ax=gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off
    set(gca,'color','none')
    print(['hrf_canon2d_ds' num2str(nd) '_' seed_list(r).name(1:end-7)], '-dpng','-r800')
    
    close all
    end
end





end %function