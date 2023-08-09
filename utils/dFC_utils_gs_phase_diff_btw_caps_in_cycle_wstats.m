function [gs_phase_diff_btw_caps_in_cycle] = dFC_utils_gs_phase_diff_btw_caps_in_cycle(phase, cfc, cap_ind, thr)

%
% This function finds and builds a distribution of global signal phase-difference at
% each CAP-pair's occurrence in two ways: conditioning the CAP's normalized cfc
% to be above a given threshold, and without the restriction. It also
% guarantees that a CAP is in a GS cycle and its pair is also inside it, or
% the preceeding or subsequent cycle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   phase: Nsubs cell with the GS phase timecourse.

%   cfc : normalized cfc timecourses (Nsubs x Ncaps).
%   cap_ind: timecourse of cap_indexes for each subject.
%  
% OUTPUTS
%   gs_phase_diff_btw_caps_in_cycle: structure with the distributions of GS phase at each
%   CAP, and circular stats.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [gs_phase_diff_btw_caps_in_cycle] = caps_analysis_gs_phase_diff_btw_caps_in_cycle(phase, cfc, cap_ind, thr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17

[Nsubs, k] = size(cfc);
for sub = 1:Nsubs
    Nobs(sub) = length(cap_ind{sub});
end



%  Detect GS cycles and their intervals (0,2pi).
for sub = 1:Nsubs
    int = 0;
    for t = 1:Nobs(sub)-1
        if sign(phase{sub}(t))==-1 && sign(phase{sub}(t+1))==1
            int = int+1;
            t_int{sub}(int) = t+1;
           
        end
    end
   
    for int = 1:length(t_int{sub})-1
            t_intr{sub}{int} = t_int{sub}(int):1:t_int{sub}(int+1);
    end   
    
end
%  Retain cycles within a range
t_min = 75;
t_max = 100;
t_intr_ths = t_intr;
for sub = 1:Nsubs
    for int = 1:length(t_intr{sub})
        if length(t_intr{sub}{int}) > t_max
            t_intr_ths{sub}{int} = [];
        elseif length(t_intr{sub}{int}) < t_min
            t_intr_ths{sub}{int} = [];
        end
        
    end
   t_intr_ths{sub}= t_intr_ths{sub}(~cellfun('isempty',t_intr{sub}));
end


% 3. Make, ofr each subject, a vector that informs if a CAP is in a determined cycle.
%Collect GS phases inside these cycles.
for c = 1:k
     gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c} = cell(Nsubs);
    for sub = 1:Nsubs
        tt=0;
        tt1=0;
        CAP_in_GSint{sub}=zeros(length(t_intr{sub}),k);
        for int = 1:length(t_intr{sub})
            for t = 1:Nobs(sub)
                if (min(t_intr{sub}{int}) < t) && (t <= max(t_intr{sub}{int})) && (cap_ind{sub}(t) == c)
                    tt=tt+1;
                    CAP_in_GSint{sub}(int,c) = 1;
                    if cfc{sub,c}(t) > thr
                        tt1=tt1+1;
                        CAP_in_GSint_ths{sub}(int,c) = 1;
                    end
                end
            end
        end
    end
end


% For each subject, and interval: if CAP1 is there, collect the time
% indexes of its presence. Then if CAP2 is present in the previous, current
% or successive cycle, collect the time indexes of its presence. Then
% compute the GS phase difference at those time indexes.
 for c1 = 1:k
    for c2 = 1:k
        
            for sub = 1:Nsubs
            tt = 0;
            tt1= 0;    
                for int = 1:length(t_intr{sub})
                    % enter each interval and ask if C1 and C2 are there
                    if CAP_in_GSint{sub}(int,c1)==1
                        %enter the cycle and mark its range
                        t_range = t_intr{sub}{int};
                        if int>1
                            t_range_low = t_intr{sub}{int-1};
                        else
                            t_range_low = [];
                        end
                        if int < length(t_intr{sub})
                            t_range_up = t_intr{sub}{int+1};
                        else
                            t_range_up = [];
                        end
                        % capture all of the occurrences of CAP1 in the
                        % interval (time indexes).
                        c1_tt = t_range(cap_ind{sub}(t_range)==c1);
                        
                        % capture all of the occurrences of CAP2 in the
                        % current, previous and posterior cycles.
                        c2_tt = [t_range_low(cap_ind{sub}(t_range_low)==c2), t_range(cap_ind{sub}(t_range)==c2), t_range_up(cap_ind{sub}(t_range_up)==c2)];
                       
                        % now compute the phase differences between the
                        % time indexes of C1 and C2.
                        if isempty(c1_tt)==0 && isempty(c2_tt)==0
                            for cc1 = 1:length(c1_tt)
                                for cc2 = 1:length(c2_tt)
                                tt=tt+1;
                                gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}{sub}(tt) = (angleDiff(phase{sub}(c2_tt(cc2)),phase{sub}(c1_tt(cc1))));
                                
                                    if cfc{sub,c1}(c1_tt(cc1)) > thr && cfc{sub,c2}(c2_tt(cc2)) > thr
                                        tt1=tt1+1;
                                        gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}{sub}(tt1) = (angleDiff(phase{sub}(c2_tt(cc2)),phase{sub}(c1_tt(cc1))));
                                    end
                                end
                            end
                        end                                                                                     
                       
                    end
                    
                end

            end
            
        
    end
 end


% compute circ stats.
% 
% circ stats for GS phase diff between CAP occurrences.
for c1 = 1:k
    for c2 = 1:k
        [gs_phase_diff_btw_caps_in_cycle.stats_thr.cmean(c1,c2),gs_phase_diff_btw_caps_in_cycle.stats_thr.cmean_CI_UP(c1,c2), gs_phase_diff_btw_caps_in_cycle.stats_thr.cmean_CI_LOW(c1,c2)]  = circ_mean(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats_thr.cvar(c1,c2) = circ_var(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats_thr.kappa(c1,c2) = circ_kappa(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats_thr.rvector(c1,c2) = circ_r(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        [gs_phase_diff_btw_caps_in_cycle.stats_thr.cstd(c1,c2),gs_phase_diff_btw_caps_in_cycle.stats_thr.cstd0(c1,c2)] = circ_std(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));

        [gs_phase_diff_btw_caps_in_cycle.stats.cmean(c1,c2),gs_phase_diff_btw_caps_in_cycle.stats.cmean_CI_UP(c1,c2), gs_phase_diff_btw_caps_in_cycle.stats.cmean_CI_LOW(c1,c2)]  = circ_mean(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats.cvar(c1,c2) = circ_var(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats.kappa(c1,c2) = circ_kappa(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats.rvector(c1,c2) = circ_r(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        [gs_phase_diff_btw_caps_in_cycle.stats.cstd(c1,c2),gs_phase_diff_btw_caps_in_cycle.stats.cstd0(c1,c2)] = circ_std(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
       
        % test for circular uniformity
        [gs_phase_diff_btw_caps_in_cycle.stats_thr.rtest_p(c1,c2), gs_phase_diff_btw_caps_in_cycle.stats_thr.rtest_z(c1,c2)] = circ_rtest(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        [gs_phase_diff_btw_caps_in_cycle.stats_thr.otest_p(c1,c2), gs_phase_diff_btw_caps_in_cycle.stats_thr.otest_m(c1,c2)] = circ_otest(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats_thr.raotest_p(c1,c2) = circ_raotest(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint_thr{c1,c2}));

         [gs_phase_diff_btw_caps_in_cycle.stats.rtest_p(c1,c2), gs_phase_diff_btw_caps_in_cycle.stats.rtest_z(c1,c2)] = circ_rtest(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        [gs_phase_diff_btw_caps_in_cycle.stats.otest_p(c1,c2), gs_phase_diff_btw_caps_in_cycle.stats.otest_m(c1,c2)] = circ_otest(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
        gs_phase_diff_btw_caps_in_cycle.stats.raotest_p(c1,c2) = circ_raotest(spm_vec(gs_phase_diff_btw_caps_in_cycle.GS_phase_diff_bw_CAPS_in_GSint{c1,c2}));
     
    end
end

%end function
end