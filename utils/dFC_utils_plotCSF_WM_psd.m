function dFC_utils_plotCSF_WM_psd


clc, clear

main = '/media/DATA/dgutierrez/human_data/caps_hnu_10session_GM_concat_3mm/analysis_CSF_WM/';
maindata = '/media/DATA/dgutierrez/human_data/HNU/test/csf_wm_3mm/';
cd(maindata)
TR(1:10)=2;
NW=4;
for nd = 1:10
    for sub = 1:30
        csf{nd}{sub} = load(['hnu_s' num2str(sub) '_ss' num2str(nd) '_vs.txt']);
        wms{nd}{sub} = load(['hnu_s' num2str(sub) '_ss' num2str(nd) '_vs.txt']);
        
        csf_psd{nd} = dFC_utils_caps_psd(csf{nd}, TR(nd));
        wms_psd{nd} = dFC_utils_caps_psd(wms{nd}, TR(nd));
        
        csf_psdMT{nd} = dFC_utils_caps_psdMT(csf{nd}, TR(nd),NW);
        wms_psdMT{nd} = dFC_utils_caps_psdMT(wms{nd}, TR(nd),NW);
        
    end
end

cd(main)
for nd = 1:10
    nsubs(nd) = length(csf{nd});
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
        fig.Position(4) = 6; fig.Position(1) = 5; fig.Position(2) = 5;   
        shadedErrorBar(csf_psdMT{nd}.min_freq_vec(10:end), csf_psdMT{nd}.psd_ds_mean{1}(10:end),csf_psdMT{nd}.psd_ds_std{1}(10:end)/sqrt(nsubs(nd)),'b')

        legend({['csf']},'FontSize', 10)
        xlabel('Frequency (Hz)', 'FontSize', 10)
        ylabel(['PSD (mean +/-SEM)'], 'FontSize', 10)
        box off
        fig.PaperPositionMode = 'auto';
        print(['psdMT_csf_ds' num2str(nd)], '-dpng','-r600')
        close all
        
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 8;
        fig.Position(4) = 6; fig.Position(1) = 5; fig.Position(2) = 5;   
        shadedErrorBar(wms_psdMT{nd}.min_freq_vec(10:end), wms_psdMT{nd}.psd_ds_mean{1}(10:end),wms_psdMT{nd}.psd_ds_std{1}(10:end)/sqrt(nsubs(nd)),'r')

        legend({['wm']},'FontSize', 10)
        xlabel('Frequency (Hz)', 'FontSize', 10)
        ylabel(['PSD (mean +/-SEM)'], 'FontSize', 10)
        box off
        fig.PaperPositionMode = 'auto';
        print(['psdMT_wm_ds' num2str(nd)], '-dpng','-r600')
        close all
end





end % function