function [trans_prob] = dFC_utils_caps_analysis_trans_prob_peaks(peaks_sub)

% This function computes the transition probability matrix from the
% CFC peaks at each CAP, condensing them into a markov_chain for each subject. It is done in two ways:
% first, taking into account self-transitions; then, without taking them
% into account.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   peaks_sub:  cell (Nsubjucts) with the peaks of CAP CFC
%   (see caps_analysis_CFC_peaks.m).

% OUTPUTS
%   trans_prob_peaks: structure witht the transition prob. matrices for each
%   subject and mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE EXAMPLE:
% [trans_prob_peaks] = caps_analysis_trans_prob_peaks(peaks_sub)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2016, V5 17-11-17

k = max(spm_vec(peaks_sub{:}));
Nsubs = length(peaks_sub);
for sub = 1:Nsubs
    Npeaks(sub) = length(peaks_sub{sub});
end

for sub = 1:Nsubs
    for c1 = 1:k
        for c2 = 1:k
            count_wa = 0;
            count_woa=0;
            % t_prob counts.
            from = peaks_sub{sub}(:);
            to   = peaks_sub{sub}(:);
            for t = 2:Npeaks(sub)
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

for sub = 1:Nsubs
    twa(sub,:,:) = trans_prob.w_auto{sub};
    twoa(sub,:,:) = trans_prob.wo_auto{sub};
end
for c1 = 1:k
    for c2 = 1:k
        trans_prob.wo_auto_mean(c1,c2) = mean(twoa(:,c1,c2));
        trans_prob.wo_auto_std(c1,c2) = std(twoa(:,c1,c2));
        trans_prob.w_auto_mean(c1,c2) = mean(twa(:,c1,c2));
        trans_prob.w_auto_std(c1,c2) = std(twa(:,c1,c2));
    end
end

%end function
end