function dFC_GS_psd_from_data3

% this function takes .nii files from a directory and computes the GS power
% spectra within a given mask.

clc,clear
% define the data directories (if more than 1)
main = '/media/DATA/dgutierrez/human_data/HNU/test/registration_3mm/';
ndatasets = 10;
ss = {'01'};
for s = 1:ndatasets
    datadir{s} = ['/media/DATA/dgutierrez/human_data/HNU/test/registration_3mm/'];
    suffix{s} = ['*ss' num2str(s) '*registered_bpf*'];
    data_name{s} = ['reg_bpf' num2str(s)];
    TR{s} = 2;
end
%define the rois to be extracted and directories
roi_dir = '/media/DATA/dgutierrez/human_data/HNU/test/registration_3mm/masks/';
nrois = 3;
roi_list_suff = {'gray_bin', 'white_bin', 'csf_bin'}; %unique suffixes
roi_name = {'gs', 'wms', 'csf'}; %unique suffixes
cd(roi_dir)
for r = 1:nrois
    rf = (dir(['*' roi_list_suff{r} '*']));
    rfile{r} = [rf.folder,'/', rf.name];
end
    
%extract data within rois
for nd = 1:ndatasets
    cd(datadir{nd}),
    list = dir(suffix{nd});
    Nsubs = length(list);
    for sub = 1:Nsubs
        sub_list{sub,1} = list(sub).name;
    end
    % Load each subject and convert dta to matrix form, and compute the gs
    % signal, and psd.
    for sub=1:Nsubs
        % Go to the subjects preprocessing dir
        for r = 1:nrois
           ts(sub,r,:) = rex(list(sub).name, rfile{r}, 'summary_measure', 'mean');
        end
            
    end
    NW = 4;
    gs_psdMT{nd} = dFC_utils_caps_psdMT(gs{nd}, TR{nd},NW);
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
    print([data_name{nd} '_gs_psdMT_ds' num2str(nd)], '-dpng','-r600')
    close all
end

end % function