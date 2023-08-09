function [motion2CAP] = dFC_utils_motion2cap(motion_info, frame_ind_sub, rm_subjects)

% This function takes the motion information from each subject and the frames associated to each CAP in that subject, and checks to which CAP do high-motion frames belong to.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   motion_info  : structure with the motion information from each the
%   dataset (see motion_check.m)
%   cap_occ      : structure with the information on CAP occurrence for
%   each subject for a specific number of clusters (see
%   caps_analysis_sub_occ.m)
% OUTPUTS
%   motion2CAP   : structure with the information of motion associated to each CAP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
%[motion2CAP] = caps_analysis_motion2CAP(motion_info, cap_occ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 15-11-17

% detect the amount of subjects and observations, and amount of clusters.
Nsubs = length(frame_ind_sub);
for sub = 1:Nsubs
    Nobs(sub) = length(frame_ind_sub{sub});
end
k = max(spm_vec(frame_ind_sub));
if isempty(rm_subjects) == 0
    motion_info.FD_flag(rm_subjects) = [];
    motion_info.FD_flag_proportion(rm_subjects) = [];
    motion_info.FD(rm_subjects) = [];
end

%% 1.associate motion volumes to its respective cap.
motion2CAP = [];
motion2CAP.count = cell(Nsubs);
motion2CAP.FD_flag_at_cap = cell(Nsubs);
motion2CAP.FD_at_cap = cell(Nsubs);
for sub = 1:Nsubs
    motion2CAP.count{sub}    = zeros(Nobs(sub),k);
    motion2CAP.FD_flag_at_cap{sub}    = zeros(Nobs(sub),k);
    motion2CAP.FD_at_cap{sub}    = zeros(Nobs(sub),k);
    fr = 0;
    for t = 1:Nobs(sub)
        if motion_info.FD_flag{sub,1}(t) == 1
            fr = fr+1;
            for c = 1:k
                if frame_ind_sub{sub}(fr) == c
                    motion2CAP.count{sub}(t,c) = 1;
                    motion2CAP.FD_flag_at_cap{sub}(t,c) = motion_info.FD{sub}(t);                    
                end                
            end            
        end
        for c = 1:k
                if frame_ind_sub{sub}(t) == c
                    motion2CAP.FD_at_cap{sub}(t,c) = motion_info.FD{sub}(t);                    
                end                
        end
    end   
    for c = 1:k
        motion2CAP.proportion_allTS(sub,c) = sum(motion2CAP.count{sub}(:,c)/Nobs(sub));
        motion2CAP.cap_proportion(sub,c) = sum(motion2CAP.count{sub}(:,c))/sum(frame_ind_sub{sub}==c);
        motion2CAP.FD_flag_at_cap_mean(sub,c)=mean(motion2CAP.FD_flag_at_cap{sub}(:,c));
        motion2CAP.FD_at_cap_mean(sub,c)=mean(motion2CAP.FD_at_cap{sub}(:,c));
        motion2CAP.count_all(sub,c)  = sum(motion2CAP.count{sub}(:,c));
    end   
end
%end function
end
