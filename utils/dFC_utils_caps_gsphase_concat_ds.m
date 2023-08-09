function dFC_utils_caps_gsphase_concat_ds

% this is an auxiliary function to re-assess GS phase distributions after
% concatenating samples from a selected group of datasets

clc, clear
%% 1. EDIT here
% load and the datasets to be merged and results paths
k=8;
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_parcel/cap_maps_2mm_k8/caps_dynamics_k8_thr1_06-1hz';
cd(main)
load('caps_dyn.mat')

ndatasets = 3; % define the amount of datasets to compare
ds_merge = [2,3];
sub_list_ds = 20:39;
% was the data scrubbed?
scrub =1;
%Load the results from each datasets for the best run
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_parcel/cap_maps_2mm_k8/cap_results_k8')
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_parcel/cap_maps_2mm_k8/inputs')
rm_subjects = inputs.rm_subjects;

%% 2. merge the datasets' main variables
nsubsT = 0;
for nd = 1:length(ds_merge)
    nsubs(nd) = inputs.Nsubs{ds_merge(nd)};
    nsubsT = nsubsT+nsubs(nd);
end

flow = caps_dyn.filter.flow;
fup = caps_dyn.filter.fup;
cfc_thr = caps_dyn.cfc.cfc_threshold;

ss = 0;
for i = 1:length(ds_merge)
    for sub = 1:nsubs(i)
        ss = ss+1;
        gs{ss,1} = caps_dyn.gs.gs{ds_merge(i)}{sub};
        gs_norm{ss,1} = caps_dyn.gs.gs_norm{ds_merge(i)}{sub};
        gs_filt{ss,1} = caps_dyn.gs.gs_filt{ds_merge(i)}{sub};
        gs_phase{ss,1} = caps_dyn.gs.gs_phase{ds_merge(i)}.phase{sub};
        for c = 1:k
            cfc{ss,c} = caps_dyn.cfc.cfc{ds_merge(i)}.cfc{sub,c};
            cfc_norm{ss,c} = caps_dyn.cfc.cfc{ds_merge(i)}.cfc_norm{sub,c};
            cfc_filt{ss,c} = caps_dyn.cfc.cfc_filt{ds_merge(i)}{sub,c};
        end
    end
end

%% 3. recompute GS-phase dsitributions at each CAP.

[gs_phase_at_cap] = dFC_utils_caps_gs_phase_at_cap(gs_phase, cfc_filt, results_map.frame_ind_sub(sub_list_ds), cfc_thr);
[gs_phase_at_cap_wtn_cycles] = dFC_utils_gs_phase_at_cap_in_cycle(gs_phase, cfc_filt,results_map.frame_ind_sub(sub_list_ds), cfc_thr);

%plot results with and without thresholding CFC


for c = 1:k
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 5;
    fig.Position(4) = 5; fig.Position(1) = 20; fig.Position(2) = 20;
    polarhistogram(spm_vec(gs_phase_at_cap_wtn_cycles.gs_phase_at_CAP_thr{c}),20,'normalization','probability')
    pax = gca; pax.ThetaTick = [0 90 180 270] ; pax.FontSize = 12; pax.RTick = [0.05 0.1 0.15];
    pax.ThetaAxisUnits = 'radians'; pax.ThetaColor = 'b'; pax.RColor = 'k';
    pax.LineWidth = 2;
    title(['GS phase at CAP ' num2str(c)], 'FontSize',12)
    set(gcf,'color','none')
    set(gca,'color','none')
    fig.PaperPositionMode = 'auto';
    print(['gs_phase_at_dsNC_cap' num2str(c)], '-dpng','-r600')
    close all
    label{c} = ['C' num2str(c)];
    
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 5;
    fig.Position(4) = 5; fig.Position(1) = 20; fig.Position(2) = 20;
    polarhistogram(spm_vec(gs_phase_at_cap_wtn_cycles.gs_phase_at_CAP{c}),20,'normalization','probability')
    pax = gca; pax.ThetaTick = [0 90 180 270] ; pax.FontSize = 12; pax.RTick = [0.05 0.1 0.15];
    pax.ThetaAxisUnits = 'radians'; pax.ThetaColor = 'b'; pax.RColor = 'k';
    pax.LineWidth = 2;
    title(['GS phase at CAP ' num2str(c)], 'FontSize',12)
    set(gcf,'color','none')
    set(gca,'color','none')
    fig.PaperPositionMode = 'auto';
    print(['gs_phase_at_dsNC_cap' num2str(c) '_nothr'], '-dpng','-r600')
    close all
    label{c} = ['C' num2str(c)];
end






end % function
