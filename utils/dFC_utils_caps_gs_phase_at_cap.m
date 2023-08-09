function [gs_phase_at_cap] = dFC_utils_caps_gs_phase_at_cap(phase, cfc, cap_ind, thr)

% This function finds and builds a distribution of global signal phases at
% each CAP's occurrence in two ways: conditioning the CAP's normalized cfc
% to be above a given threshold, and without the restriction.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   phase: Nsubs cell with the GS phase timecourse.
%   cfc : normalized cfc timecourses (Nsubs x Ncaps).
%   cap_ind: timecourse of cap_indexes for each subject.

% OUTPUTS
%   gs_phase_at_cap: structure with the distributions of GS phase at each
%   CAP, and circular stats.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [gs_phase_at_cap] = caps_analysis_gs_phase_at_cap(phase, cfc, cap_ind, thr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17


[Nsubs, k] = size(cfc);
for sub = 1:Nsubs
    Nobs(sub) = length(cap_ind{sub});
end

% 1. collect GS phases at the appearance of each cap (with and without
% thresholding occurrences to be above a cfc threshold).

for c = 1:k
    for sub = 1:Nsubs
        tt=0;
        tt2=0;
        for t = 1:Nobs(sub)
            if cap_ind{sub}(t) == c
                tt=tt+1;
                gs_phase_at_cap.gs_phase_at_CAP{c}{sub}(tt) = phase{sub}(t);
                
                if cfc{sub,c}(t) > thr
                    tt2=    tt2+1;
                    gs_phase_at_cap.gs_phase_at_CAP_thr{c}{sub}(tt2) = phase{sub}(t);
                    
                end
                
            end
        end
    end
end


% compute circ stats.
for c = 1:k
    [gs_phase_at_cap.stats.cmean(c), gs_phase_at_cap.stats.cmean_CI_UP(c), gs_phase_at_cap.stats.cmean_CI_LOW(c)]  = circ_mean(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));
    gs_phase_at_cap.stats.cvar(c) = circ_var(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));
    gs_phase_at_cap.stats.kappa(c) = circ_kappa(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));
    gs_phase_at_cap.stats.rvector(c) = circ_r(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));
    gs_phase_at_cap.stats.cstd(c) = circ_std(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));
    [gs_phase_at_cap.stats.cstd(c),gs_phase_at_cap.stats.cstd0(c)] = circ_std(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));

    
   [gs_phase_at_cap.stats_thr.cmean(c), gs_phase_at_cap.stats_thr.cmean_CI_UP(c), gs_phase_at_cap.stats_thr.cmean_CI_LOW(c)]  = circ_mean(spm_vec(gs_phase_at_cap.gs_phase_at_CAP_thr{c}));
    gs_phase_at_cap.stats_thr.cvar(c) = circ_var(spm_vec(gs_phase_at_cap.gs_phase_at_CAP_thr{c}));
    gs_phase_at_cap.stats_thr.kappa(c) = circ_kappa(spm_vec(gs_phase_at_cap.gs_phase_at_CAP_thr{c}));
    gs_phase_at_cap.stats_thr.rvector(c) = circ_r(spm_vec(gs_phase_at_cap.gs_phase_at_CAP_thr{c}));
    [gs_phase_at_cap.stats_thr.cstd(c),gs_phase_at_cap.stats_thr.cstd0(c)] = circ_std(spm_vec(gs_phase_at_cap.gs_phase_at_CAP_thr{c}));


    %test for circular uniformity (Raleigh test)
    [gs_phase_at_cap.stats.rtest_p(c), gs_phase_at_cap.stats.rtest_z(c)] = circ_rtest(spm_vec(gs_phase_at_cap.gs_phase_at_CAP{c}));
    [gs_phase_at_cap.stats_thr.rtest_p(c), gs_phase_at_cap.stats_thr.rtest_z(c)] = circ_rtest(spm_vec(gs_phase_at_cap.gs_phase_at_CAP_thr{c}));
end

%end function
end