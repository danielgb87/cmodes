function dFC_GS_psd_from_data

% this function takes .nii files from a directory and computes the GS power
% spectra within a given mask.

clc,clear

main = '/media/DATA/dgutierrez/human_data/HNU_fullPP_Ting/smoothing/sess_01/';
ndatasets = 1;
ss = {'01'};
for s = 1:ndatasets
    datadir{s} = ['/media/DATA/dgutierrez/human_data/HNU_fullPP_Ting/smoothing/sess_01/'];
    suffix{s} = ['*ss' num2str(s) '*_fullPP*'];
    data_name{s} = ['ds1_brainfullPPTing'];
    TR{s} = 2;
end

mask_file = '/media/DATA/dgutierrez/human_data/masks_templates_parcels/fsl_masks/MNI152_T1_3mm_brain_mask.nii.gz';
% Check indexes with non-zero values
mask = single(spm_read_vols(spm_vol(mask_file)));
[mask_ind(1,:), mask_ind(2,:), mask_ind(3,:)] = ind2sub(size(mask),find(mask ==1));
mask_info.file = mask_file;
mask_info.mask = mask;
mask_info.mask_index = mask_ind;

%extract data within mask
for nd = 1:ndatasets
    cd(datadir{nd}),
    list = dir(suffix{nd});
    Nsubs = length(list);
    for sub = 1:Nsubs
        sub_list{nd}{sub,1} = list(sub).name;
    end
    % Load each subject and convert dta to matrix form, and compute the gs
    % signal, and psd.
    for sub=1:Nsubs
        % Go to the subjects preprocessing dir   
            temp = spm_vol(sub_list{nd}{sub,1});
            parfor t = 1:length(temp)
                data1(:,t) = spm_vec(spm_read_vols(temp(t)));
            end
            data{nd}{sub,1} = data1(spm_vec(mask)==1,:)';
            gs{nd}{sub,1} = zscore(mean(data{nd}{sub},2));
            clear data1 temp
    end
    NW = 4;
    gs_psdMT{nd} = dFC_utils_caps_psdMT(gs{nd}, TR{nd},NW);
    disp(['don ds' num2str(nd)])
end

%plot PSD
cd(main)
for nd = 1:ndatasets
    nsubs(nd) = length(gs{nd});
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
    fig.Position(4) = 6; fig.Position(1) = 5; fig.Position(2) = 5;
    shadedErrorBar(gs_psdMT{nd}.min_freq_vec, gs_psdMT{nd}.psd_ds_mean{1},gs_psdMT{nd}.psd_ds_std{1}/sqrt(nsubs(nd)),'b')
    [M,ii] =max(gs_psdMT{nd}.psd_ds_mean{1}); fp = gs_psdMT{nd}.min_freq_vec(ii);
    ymax = max(max(spm_vec(gs_psdMT{nd}.psd_ds_mean{1})))+1.2*max(max(spm_vec(gs_psdMT{nd}.psd_ds_std{1})))/sqrt(nsubs(nd));
    ylim([0,round(ymax,1)]), xlim([0 0.15])
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
    print(['aaGS_' data_name{nd} '_psd_ds' num2str(nd)], '-dpng','-r600')
    close all
end

end % function