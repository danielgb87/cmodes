function  [cap_cfc_analysis, btw_cap_corr] = dFC_utils_caps_cfc(data_chosen, templates)


% This function takes K fmri maps (templates) with N in-brain voxels and
% first, computes the spatial correltion between them, then using a dataset
% with M cells (subjects) with a T-time x N-voxels time-frames, and builds
% timecourses of template to frame correlations for each template and each
% subject .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   data_chosen : cell structure (one cell per subject) containing the data
%   matrix (time x voxels).
%   templates   : matrix (K-templates x N-voxels) conteining the template
%   maps.

% OUTPUTS
%   btw_cap_corr: spatial correlation between templates.
%   cap_CFC_analysis: structure with the template to frame timecourses for
%   each subject and template.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%[cap_CFC_analysis, btw_cap_corr] = caps_analysis_CFC(data_chosen,
%results);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 14-11-17

% Detect amount of templates, subjects, and frames per subject. Also detect the
% amount of clusters  (k) in the analysis.

k = size(templates,1);
Nsubs = length(data_chosen);
for sub = 1:Nsubs
    Nobs(sub) = size(data_chosen{sub},1);
end
% compute between CAP similarity matrix.
for c1 = 1:k
    for c2 = 1:k
        btw_cap_corr(c1,c2) = corr(spm_vec(templates(c1,:)),spm_vec(templates(c2,:)));
    end
end

% now compute CAP-to-frame spatial correlation timecourses for each
% subject.
for sub = 1:Nsubs
    for c = 1:k
        for t = 1:Nobs(sub)
            cap_cfc_analysis.cfc{sub,c}(t) = corr(spm_vec(templates(c,:)), spm_vec(data_chosen{sub}(t,:)));
        end        
        % normalized timecourses
        cap_cfc_analysis.cfc_norm{sub,c} = zscore(cap_cfc_analysis.cfc{sub,c}(:));
    end    
end

% end function.
end

