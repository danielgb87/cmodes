function [gs_phase_at_cap_in_cycle] = dFC_utils_gs_phase_at_cap_in_cycle(phase, cfc, cap_ind, thr)

% This function finds and builds a distribution of global signal phases at
% each CAP's occurrence in two ways: conditioning the CAP's normalized cfc
% to be above a given threshold, and without the restriction. Here however,
% the Gs is partitioned in cycles and the occurrence of CAPs is evaluated
% in each cycle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   phase: Nsubs cell with the GS phase timecourse.

%   phase: cell with the phase timecourses for each subject.
%   cfc : normalized cfc timecourses (Nsubs x Ncaps).
%   cap_ind: timecourse of cap_indexes for each subject.
%  
% OUTPUTS
%   gs_phase_at_cap_in_cycle: structure with the distributions of GS phase at each
%   CAP, and circular stats.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [gs_phase_at_cap_in_cycle] = caps_analysis_gs_phase_at_cap_in_cycle(phase, cfc, cap_ind, thr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17


[Nsubs, k] = size(cfc);
for sub = 1:Nsubs
    Nobs(sub) = length(cap_ind{sub});
end



% 2. Detect GS cycles and their intervals (0,2pi).
for sub = 1:Nsubs
    int = 0;
    for t = 1:Nobs(sub)-1
        if sign(phase{sub}(t))==-1 && sign(phase{sub}(t+1))==1
%         if cfc.phase_analysis.GS_phase(sub,t)<pi/2 && cfc.phase_analysis.GS_phase(sub,t+1)>pi/2
            int = int+1;
            t_int{sub}(int) = t+1;
           
        end
    end
   
    for int = 1:length(t_int{sub})-1
            t_intr{sub}{int} = t_int{sub}(int):1:t_int{sub}(int+1);
    end   
    
end
% 2.1 Retain cycles within a range
t_min = 30;
t_max = 70;
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
        for int = 1:length(t_intr{sub})
            for t = 1:Nobs(sub)
                if (min(t_intr{sub}{int}) < t) && (t <= max(t_intr{sub}{int})) && (cap_ind{sub}(t) == c)
                    tt=tt+1;
                    gs_phase_at_cap_in_cycle.CAP_in_GSint{sub}(int,c) = 1;
                    gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}{sub}(tt) = phase{sub}(t);
                    if cfc{sub,c}(t) > thr
                        tt1=tt1+1;
                        gs_phase_at_cap_in_cycle.CAP_in_GSint_ths{sub}(int,c) = 1;
                        gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}{sub}(tt1) = phase{sub}(t);
                    end
                end
            end
        end
    end
end



% compute circ stats.
for c = 1:k
    if isempty(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}))==0
        [gs_phase_at_cap_in_cycle.stats.cmean(c), gs_phase_at_cap_in_cycle.stats.cmean_CI_UP(c), gs_phase_at_cap_in_cycle.stats.cmean_CI_LOW(c)]  = circ_mean(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));
        gs_phase_at_cap_in_cycle.stats.cvar(c) = circ_var(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));
        gs_phase_at_cap_in_cycle.stats.kappa(c) = circ_kappa(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));
        gs_phase_at_cap_in_cycle.stats.rvector(c) = circ_r(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));
        gs_phase_at_cap_in_cycle.stats.cstd(c) = circ_std(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));
        [gs_phase_at_cap_in_cycle.stats.cstd(c),gs_phase_at_cap_in_cycle.stats.cstd0(c)] = circ_std(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));
       %test for circular uniformity (Raleigh test)
        [gs_phase_at_cap_in_cycle.stats.rtest_p(c), gs_phase_at_cap_in_cycle.stats.rtest_z(c)] = circ_rtest(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP{c}));

    end
    
    if isempty(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}))==0
    
   [gs_phase_at_cap_in_cycle.stats_thr.cmean(c), gs_phase_at_cap_in_cycle.stats_thr.cmean_CI_UP(c), gs_phase_at_cap_in_cycle.stats_thr.cmean_CI_LOW(c)]  = circ_mean(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}));
    gs_phase_at_cap_in_cycle.stats_thr.cvar(c) = circ_var(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}));
    gs_phase_at_cap_in_cycle.stats_thr.kappa(c) = circ_kappa(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}));
    gs_phase_at_cap_in_cycle.stats_thr.rvector(c) = circ_r(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}));
    [gs_phase_at_cap_in_cycle.stats_thr.rtest_p(c), gs_phase_at_cap_in_cycle.stats_thr.rtest_z(c)] = circ_rtest(spm_vec(gs_phase_at_cap_in_cycle.gs_phase_at_CAP_thr{c}));

    end
end

%end function
end