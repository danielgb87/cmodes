function [trans_prob] = dFC_utils_caps_analysis_trans_prob(frame_ind)

% This function computes the transition probability matrix from the
% cap_indexes of each frame, for each subject. It is done in two ways:
% first, taking into account self-transitions; then, without taking them
% into account.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   frames_ind:  cell (N subjects) with the timecourse of CAP occurrence
%   (see caps_analysis_sub_occ.m).

% OUTPUTS
%   trans_prob: structure witht the transition prob. matrices for each
%   subject and mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17

k = max(spm_vec(frame_ind{:}));
Nsubs = length(frame_ind);
for sub = 1:Nsubs
    Nobs(sub) = length(frame_ind{sub});
end

for sub = 1:Nsubs
    for c1 = 1:k
        for c2 = 1:k
            count_wa = 0;
            count_woa = 0;
            % t_prob counts.
            from = frame_ind{sub}(:);
            to   = frame_ind{sub}(:);
            for t = 2:Nobs(sub)
                if from(t-1) ==c1 && to(t) == c2
                    count_wa = count_wa+1;
                end
                if c1 ~= c2
                    if from(t-1) ==c1 && to(t) == c2
                        count_woa = count_woa+1;
                    end
                end
            end
            
            trans_prob.w_auto{sub}(c2,c1) = count_wa;
            trans_prob.wo_auto{sub}(c2,c1) = count_woa;
            
        end
    end
    
    %normalize by columns
    for c = 1:k
        trans_prob.w_auto{sub}(:,c) = trans_prob.w_auto{sub}(:,c)/sum(trans_prob.w_auto{sub}(:,c));
        trans_prob.wo_auto{sub}(:,c) = trans_prob.wo_auto{sub}(:,c)/sum(trans_prob.wo_auto{sub}(:,c));
    end
end



%Reorganize for descriptive stats

for sub = 1:size(frame_ind,1)
    twa(sub,:,:) = trans_prob.w_auto{sub};
    twoa(sub,:,:) = trans_prob.wo_auto{sub};
end
for c1 = 1:k
    for c2 = 1:k
        trans_prob.w_auto_mean(c1,c2) = mean(twa(:,c1,c2));
        trans_prob.wo_auto_mean(c1,c2) = mean(twoa(:,c1,c2));
        trans_prob.w_auto_std(c1,c2) = std(twa(:,c1,c2));
        trans_prob.wo_auto_std(c1,c2) = std(twoa(:,c1,c2));
    end
end

%end function
end