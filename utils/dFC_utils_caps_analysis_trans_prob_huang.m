function [count, tp] = dFC_utils_caps_analysis_trans_prob_huang(frame_ind)

% This function computes the transition count and probability matrix from the
% cap_indexes of each frame, as the probability of transitioning from A to
% B, given that the current state is A. That is, the number of transitions
% from A to B divided by the number of occurrences of A.

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

k = max(spm_vec(frame_ind));
Nobs = length(frame_ind);

count=zeros(k,k);
tp=zeros(k,k);
for c1 = 1:k
    for c2 = 1:k
        wa = 0;
        % t_prob counts.
        from = frame_ind;
        to   = frame_ind;
        for t = 2:Nobs
            if from(t-1) ==c1 && to(t) == c2
                wa = wa+1;
            end
        end        
        count(c1,c2) = wa;
        tp(c1,c2)=count(c1,c2)/(sum(frame_ind==c1));
    end
end

%end function
end