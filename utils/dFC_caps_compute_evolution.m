function dFC_caps_compute_evolution

% this function computes the CAP's evolution at the vicinity of a peak of its timecourse.
% REEMEMBER THAT THE ORDER OF THE CAPS HERE IS
% THE SAME AS THE RESULTS FILE, AND AS THEY AREMAPPED. BUT THE ORDER OF THE
% PRELIMINARY "CAPS_PROPERTIES", "CAPS_QUALITY", "CAPS_DYNAMICS" ARE DIFFERENT AND
% ARRANGED ACCOORDING TO THE CAP OCCURRENCE RATE FROM THE CONCATENATED
% ANALYSIS. 

clc,clear
%% EDIT HERE
%% 1.a load caps results, inputs, mask, and parcellation info.
k=8;
main = ['/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_ds/cap_maps_2mm_k' num2str(k)];
main_dyn = [main '/caps_dynamics_k8'];
cd(main)
ndatasets = 3; % define the amount of datasets to compare
ncaps = 2:20;    % the clusterings you want to compare... make sure all datasets and runs have them.
main_ds = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/';
dataset_name = 'concat_vw_ds';
nruns = 5;
% was the data scrubbed?
scrub =1;

%Load the results from each datasets for the best run
cd(main)
load(['cap_results_k' num2str(k)])
load('inputs')
rm_subjects = inputs.rm_subjects;

% load the data to compute VE from, and rename it
cd(main_ds)
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

cd(main_dyn)
load('caps_dyn')

TR=inputs.analysis_TR;

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

%define the threshold for peaks and the range to compute CAP evo in TRs,
%and percent of highest samples to take
cfc_thr = 1.5;
perc = 0.2;
range = 15;
%correct the FD flags if data was scrubbed
if scrub ==1
   for nd = 1:ndatasets
       for sub = 1:nsubs(nd)
        motion_info{nd}.FD_flag{sub} = zeros(size(data{nd}{sub},1),1);
       end
   end
end

%% 2. Find the CAP peaks
cd(main)
for nd = 1:ndatasets
    [peak_ind{nd}, peak_amp{nd}] = dFC_utils_caps_cfc_peaks(caps_dyn.cfc.cfc{nd}.cfc_norm, caps_dyn.cfc.cfc{nd}.cfc, motion_info{nd}, cfc_thr);   
    [evo{nd}] = dFC_utils_cap_evo_group(data{nd}, peak_ind{nd}, peak_amp{nd}, range,caps_dyn.cfc.cfc{nd}.cfc,perc);
end

%% 3. now map the CAP evolution sequence.
mkdir('cap_evo_15thr_20p')
cd('cap_evo_15thr_20p')
for nd = 1:ndatasets
    [Nstep, Nvox] = size(evo{nd}.frames_mean{1});
    [aa,bb,cc]=size(mask_info.mask);
    ind = mask_info.mask_index;
    for c = 1:k
        map = zeros(aa,bb,cc,Nstep);
        mapT = zeros(aa,bb,cc,Nstep);
        map0 = zeros(aa,bb,cc);
        for t = 1:Nstep
            for vox = 1:Nvox
                map(ind(1,vox),ind(2,vox),ind(3,vox),t)  = evo{nd}.frames_mean{c}(t,vox);
                mapT(ind(1,vox),ind(2,vox),ind(3,vox),t)  = evo{nd}.frames_mean_T{c}(t,vox);
            end
        end
        V = spm_vol(mask_info.file);
        V.fname = [pwd '/cap_evo_mean0_ds' num2str(nd) '_c' num2str(c) '.nii'];
        V.dt=[16,0];
        spm_write_vol(V,squeeze(map(:,:,:,floor(Nstep/2+1))));
             
        t0=load_nii([pwd '/cap_evo_mean0_ds' num2str(nd) '_c' num2str(c) '.nii']);
        t1=t0; t1.img=[]; t1.img=single(map(fliplr(1:aa),:,:));
        t1.hdr.dime.dim = [4 aa bb cc Nstep 1 1 1];
        save_nii(t1, [pwd '/cap_evo_meanAssembly_ds' num2str(nd) '_c' num2str(c) '.nii'])
        gzip([pwd '/cap_evo_meanAssembly_ds' num2str(nd) '_c' num2str(c) '.nii'])
        delete ([pwd '/cap_evo_meanAssembly_ds' num2str(nd) '_c' num2str(c) '.nii'])
        
        tt0=load_nii([pwd '/cap_evo_mean0_ds' num2str(nd) '_c' num2str(c) '.nii']);
        tt1=tt0; tt1.img=[]; tt1.img=single(mapT(fliplr(1:aa),:,:));
        tt1.hdr.dime.dim = [4 aa bb cc Nstep 1 1 1];
        save_nii(tt1, [pwd '/cap_evo_meanAssembly_T_ds' num2str(nd) '_c' num2str(c) '.nii'])
        gzip([pwd '/cap_evo_meanAssembly_T_ds' num2str(nd) '_c' num2str(c) '.nii'])
        delete ([pwd '/cap_evo_meanAssembly_T_ds' num2str(nd) '_c' num2str(c) '.nii'])
        
        gzip([pwd '/cap_evo_mean0_ds' num2str(nd) '_c' num2str(c) '.nii'])
        delete ([pwd '/cap_evo_mean0_ds' num2str(nd) '_c' num2str(c) '.nii'])
    end
end


%% 4. plot the CAP evolution kernels
for nd = 1:ndatasets
    for c1 = 1:k
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 12;
        fig.Position(4) = 9; fig.Position(1) = 5; fig.Position(2) = 5;   
        shadedErrorBar(TR*evo{nd}.steps, evo{nd}.cfc_samples_mean{c1},  evo{nd}.cfc_samples_sem{c1},'-b',1)
        xlim([-TR*range TR*range])
        hAx=gca;
        set(hAx,'xtick',-TR*range:3*TR:TR*range)
        hAx.XAxisLocation = 'origin';
        hAx.YAxisLocation = 'origin';
        legend({['CAP ' num2str(c1)]},'FontSize', 10)
        xlabel('time (s)', 'FontSize', 10)
        box off,  fig.PaperPositionMode = 'auto';
        print(['cap_evo_ds' num2str(nd) '_c' num2str(c1)], '-dpng','-r600')
        close all
    end    
end


%% 5. ORGANIZE INFORMATION IN A STRUCTURE
cap_evo.evo = evo;
cap_evo.cfc_thr = cfc_thr;
cap_evo.peak_ind = peak_ind;
cap_evo.peak_amp = peak_amp;
cap_evo.range = range;
save('cap_evo','cap_evo')


end % function