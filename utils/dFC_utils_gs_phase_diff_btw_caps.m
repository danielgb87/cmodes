function [gs_phase_diff_btw_caps] = dFC_utils_gs_phase_diff_btw_caps(phase, cfc, cap_ind, thr)

%
% This function finds and builds a distribution of global signal phase-difference at
% each CAP-pair's occurrence in two ways: conditioning the CAP's normalized cfc
% to be above a given threshold, and without the restriction.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   phase: Nsubs cell with the GS phase timecourse.

%   cfc : normalized cfc timecourses (Nsubs x Ncaps).
%   cap_ind: timecourse of cap_indexes for each subject.
%  
% OUTPUTS
%   gs_phase_diff_btw_caps: structure with the distributions of GS phase at each
%   CAP, and circular stats.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   [gs_phase_diff_btw_caps] = caps_analysis_gs_phase_diff_btw_caps(phase, cfc, cap_ind, thr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17

[Nsubs, k] = size(cfc);
for sub = 1:Nsubs
    Nobs(sub) = length(cap_ind{sub});
end


% 5. Now compute GS phase differences between CAP occurrences within the
 % whole session, disregarding the cycle they pertain, but keeping with the
 % cfc threshold restriction, and without.
 
 for c1 = 1:k
     for c2 = 1:k
         
         for sub = 1:Nsubs
             tt = 0;
             tt1=0;
             for t1 = 1:Nobs(sub)-1
                 if cap_ind{sub}(t1)==c1 && cfc{sub,c1}(t1) > thr
                     for t2 = 1:Nobs(sub)
                         if cap_ind{sub}(t2)==c2 && cfc{sub,c2}(t2) > thr
                             
                             tt=tt+1;
                             gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}{sub}(tt) = (angleDiff(phase{sub}(t2),phase{sub}(t1)));
                         end
                         
                     end
                 end
             end
             
             for t1 = 1:Nobs(sub)-1
                 if cap_ind{sub}(t1)==c1
                     for t2 = 1:Nobs(sub)
                         if  cap_ind{sub}(t2)==c2
                             
                             tt1=tt1+1;
                             gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}{sub}(tt1) = (angleDiff(phase{sub}(t2),phase{sub}(t1)));
                         end
                         
                     end
                 end
             end
             
             
             
         end
     end
 end
 
% compute circ stats.

% circ stats for GS phase diff between CAP occurrences.
for c1 = 1:k
    for c2 = 1:k
        [gs_phase_diff_btw_caps.stats_thr.cmean(c1,c2),gs_phase_diff_btw_caps.stats_thr.cmean_CI_UP(c1,c2), gs_phase_diff_btw_caps.stats_thr.cmean_CI_LOW(c1,c2)]  = circ_mean(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps.stats_thr.cvar(c1,c2) = circ_var(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps.stats_thr.kappa(c1,c2) = circ_kappa(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps.stats_thr.rvector(c1,c2) = circ_r(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        [gs_phase_diff_btw_caps.stats_thr.cstd(c1,c2),gs_phase_diff_btw_caps.stats_thr.cstd0(c1,c2)] = circ_std(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));

        [gs_phase_diff_btw_caps.stats.cmean(c1,c2),gs_phase_diff_btw_caps.stats.cmean_CI_UP(c1,c2), gs_phase_diff_btw_caps.stats.cmean_CI_LOW(c1,c2)]  = circ_mean(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps.stats.cvar(c1,c2) = circ_var(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps.stats.kappa(c1,c2) = circ_kappa(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps.stats.rvector(c1,c2) = circ_r(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        [gs_phase_diff_btw_caps.stats.cstd(c1,c2),gs_phase_diff_btw_caps.stats.cstd0(c1,c2)] = circ_std(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
       
        % test for circular uniformity
        [gs_phase_diff_btw_caps.stats_thr.rtest_p(c1,c2), gs_phase_diff_btw_caps.stats_thr.rtest_z(c1,c2)] = circ_rtest(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        [gs_phase_diff_btw_caps.stats_thr.otest_p(c1,c2), gs_phase_diff_btw_caps.stats_thr.otest_m(c1,c2)] = circ_otest(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps.stats_thr.raotest_p(c1,c2) = circ_raotest(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));

         [gs_phase_diff_btw_caps.stats.rtest_p(c1,c2), gs_phase_diff_btw_caps.stats.rtest_z(c1,c2)] = circ_rtest(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        [gs_phase_diff_btw_caps.stats.otest_p(c1,c2), gs_phase_diff_btw_caps.stats.otest_m(c1,c2)] = circ_otest(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps.stats.raotest_p(c1,c2) = circ_raotest(spm_vec(gs_phase_diff_btw_caps.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
     
    end
end

%end function
end