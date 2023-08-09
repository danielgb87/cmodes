function dFC_utils_cap_evo_images(cap_evol, mask_info, prefix)


% This function takes the frames from peaks, before them and after them,
% and builds a succession of fMRI maps for each CAP. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   cap_evol:  structure containing the averaged frames in each time-step.
%   mask_info: structure containing the volume information of the brain
%   mask used thorughtout the analysis.

% OUTPUTS
% CAP evolution fMRI maps (sequence) saved as .nii files in the current
% directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
% [evo_group] = caps_analysis_cap_evo_group(data_chosen, peak_ind, peak_amp, occ_prob,range, CFC, perc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 16-11-17

k = length(cap_evol.frames_mean);
[Nstep, Nvox] = size(cap_evol.frames_mean{1});
[aa,bb,cc]=size(mask_info.mask);

ind = mask_info.mask_index;

for c = 1:k
    
    map = zeros(aa,bb,cc,Nstep);
    map0 = zeros(aa,bb,cc);
    
    for t = 1:Nstep
        
        for vox = 1:Nvox
            map(ind(1,vox),ind(2,vox),ind(3,vox),t)  = cap_evol.frames_mean{c}(t,vox);
        end
        
    end
    
    
    V = spm_vol(mask_info.file);
    V.fname = [pwd '/' prefix '_CAP_' num2str(c) '_mean_assembly_0.nii'];
    V.dt=[16,0];
    spm_write_vol(V,squeeze(map(:,:,:,floor(Nstep/2+1))));
    
    t0=load_nii([pwd '/' prefix '_CAP_' num2str(c) '_mean_assembly_0.nii']);
    t1=t0; t1.img=[]; t1.img=single(map(fliplr(1:aa),:,:));
    t1.hdr.dime.dim = [4 aa bb cc Nstep 1 1 1];
    save_nii(t1, [prefix '_CAP_' num2str(c) '_mean_assembly.nii'])
    
end

cd ..

end














