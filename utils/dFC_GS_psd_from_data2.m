function dFC_GS_psd_from_data2

% this function takes .nii files from a directory and computes the GS power
% spectra within a given mask.

clc,clear

main = '/media/DATA/dgutierrez/human_data/MSC/MSC_rsfMRI_pp2/06_filtering_smoothing/smoothing';
ndatasets = 5;
ss = {'01','02','05','08','09'};
for s = 1:ndatasets
    datadir{s} = ['/media/DATA/dgutierrez/human_data/MSC/MSC_rsfMRI_pp2/06_filtering_smoothing/smoothing/sess_' ss{s}];
    suffix{s} = ['*smoothed.nii*'];
    data_name{s} = ['pp_smoothed_ds' ss{s}];
    TR{s} = 2.2;
end

mask_file = '/media/DATA/dgutierrez/human_data/masks_templates_parcels/fsl_masks/BN_atlas_mask_3mm.nii.gz';
% Check indexes with non-zero values
for i = 1:ndatasets
    [data{i}, mask_info, sub_list{i}] = dFC_utils_prepare_dataset_matrix(datadir{i},suffix{i}, mask_file);
    nsubs(i) = length(data{i});
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
    NW = 4;
%     gs_psdMT{nd} = dFC_utils_caps_psdMT(gs_norm{nd},  TR{nd},NW);
    gs_psdMT{nd} = dFC_utils_caps_psdMT(gs_norm{nd},  TR{nd},NW);
end

% plot PSD
for nd = 1:ndatasets
    cd(datadir{nd})
    % now GS PSD
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
    fig.Position(4) = 6; fig.Position(1) = 5; fig.Position(2) = 5;
    shadedErrorBar(gs_psdMT{nd}.min_freq_vec, gs_psdMT{nd}.psd_ds_mean{1},gs_psdMT{nd}.psd_ds_std{1}/sqrt(nsubs(nd)),'b')
    [M,ii] =max(gs_psdMT{nd}.psd_ds_mean{1}); fp = gs_psdMT{nd}.min_freq_vec(ii);
    ymax = max(max(spm_vec(gs_psdMT{nd}.psd_ds_mean{1})))+1.2*max(max(spm_vec(gs_psdMT{nd}.psd_ds_std{1})))/sqrt(nsubs(nd));
    xlim([0 0.15])
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
    print(['aGS_psd_' data_name{nd}], '-dpng','-r600')
    close all
end


end % function