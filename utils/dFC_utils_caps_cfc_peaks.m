function  [peak_ind, peak_amp] = dFC_utils_caps_cfc_peaks(cfc_norm, cfc, motion_info, thr)


% This function selects for each subject, and CAP, frame indexes with the high
% quality clustering (CFC-normalized peaks with derivative zero and
% CFC_norm > thr). It assures that frames are of CFC > thr and local peaks; that they are
% not associated to high-motion; and that the CAP has the highest CFC out
% of all at that specific time-point.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   CFC      : cell matrix (Nsubjects x NCAPs) with normalized CFC timecourses
%   motion_info : structure with the motion information for the dataset.
%   thr: threshold in SD units that CFC's must be above to be selected.

% OUTPUTS
%   peak_ind  : a cell (Nsubject x NCAPs) with the highest peaks.
%   peak_amp  : a cell (Nsubject x NCAPs) with the amplitude of the highest peaks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
% [peak_ind, peak_amp] = caps_analysis_CFC_peaks(CFC,motion_info, thr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 14-11-17


mot_flags = motion_info.FD_flag;

% detect amount of subject,caps, and observations.
[Nsubs,k] = size(cfc_norm);
for sub = 1:Nsubs
    Nobs(sub) = length(cfc_norm{sub,1});
end

% first reconfigure CFC from all caps at all time points for each
% subject.
for sub = 1:Nsubs
    for c = 1:k
        cfc_temp{sub}(c,:) = cfc_norm{sub,c};
    end
end

% find local maxima (peaks)
peak_ind = cell(Nsubs,k);
peak_amp = cell(Nsubs,k);

for sub = 1:Nsubs
    for c = 1:k
        ind = 0;
        for t = 2:Nobs(sub)-1 % avoid picking frames from teh first or last timepoints.
            % check for maxima
            if cfc_norm{sub,c}(t-1) < cfc_norm{sub,c}(t) && cfc_norm{sub,c}(t+1) < cfc_norm{sub,c}(t)
                % check if CFC is above threshold.
                if cfc_norm{sub,c}(t) > thr
                    % check that time-point is not high moption.
                    if mot_flags{sub}(t) == 0
                        % check that that CAP is the highest CFC at the
                        % time-point.
                        if cfc_norm{sub,c}(t) == max(cfc_temp{sub}(:,t))
                            ind = ind+1;
                            peak_ind{sub,c}(ind) = t;
                            peak_amp{sub,c}(ind) = cfc{sub,c}(t);
                        end
                    end
                end
            end
        end
    end
    
end

%end function
end





