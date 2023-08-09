function  [evo] = dFC_utils_cap_evo_group(data_chosen, peak_ind, peak_amp, range,cfc,perc)


% This function selects for each subject, and CAP, frame indexes with the high
% quality clustering and averages them. It also averages the previous and
% subsequent frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS.
%   data_chosen      : cell matrix (Nsubjects) with data fro each subject
%   peak_ind: cell (Nsubs x Ncaps) of vectors of peak indexes
%   peak_amp:  cell (Nsubs x Ncaps) of vectors of cfc amplitudes at peaks
%   occ:  mean occurrence rate of each CAP.
%   cfc: cfc (Nsubs x Ncap) normalized timecourses.
%   thr: threshold in SD units that cfc's must be above to be selected.

% OUTPUTS
%   evo: structure with the relevant info of CAP evolution kernels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Nsubs, k] = size(cfc);
for sub = 1:Nsubs
    Nobs(sub) = length(cfc{sub,1});
end

% select the amount of samplef from the whole dataset to chose, based on
% the occ_rate.
Nobs_total = sum(Nobs);
evo.steps = -range:1:range;

% take away peaks at the first and final part of each session (depending on
% range).
for c = 1:k
    tt=0;
    for sub = 1:Nsubs
        for ind = 1:length(peak_ind{sub,c})
            if peak_ind{sub,c}(ind)+range < Nobs(sub) && peak_ind{sub,c}(ind) - range > 0
                tt=tt+1;
                
                peak_index{c}{tt} = [sub,peak_ind{sub,c}(ind)];
                peak_amps{c}(tt) = peak_amp{sub,c}(ind);
            end
        end
    end
    
    [evo.peaks_sort_amp{c},evo.peaks_sort_ind{c}] = sort(peak_amps{c},'descend');
    
    evo.peaks_ind_final{c} = peak_index{c}(evo.peaks_sort_ind{c});
    evo.nsamples(c) =ceil(length(evo.peaks_ind_final{c})*perc);
  
    % now collect the frames for each cap, their preceeding and subsequent
    % ones time-locked, and average them. Also collect the cfc amplitudes.
    for st =1:length(evo.steps)
        for fr = 1:evo.nsamples(c)
            frames{c}{st}(fr,:) = data_chosen{evo.peaks_ind_final{c}{fr}(1)}(evo.peaks_ind_final{c}{fr}(2)+evo.steps(st),:);
            evo.cfc_samples{c}{st}(fr) = cfc{evo.peaks_ind_final{c}{fr}(1),c}(evo.peaks_ind_final{c}{fr}(2)+evo.steps(st));
        end
        evo.frames_mean{c}(st,:) = mean(frames{c}{st},1);
        nvox = size(evo.frames_mean{c},2);
        temp = frames{c}{st}(:,:);
        parfor v = 1:nvox           
                [~,temp_p(v) , ~, stats] = ttest(spm_vec(temp(:,v)));
                temp_T(v) = stats.tstat;            
        end
        clear temp
        evo.frames_mean_p{c}(st,:) = temp_p;
        evo.frames_mean_T{c}(st,:) = temp_T;
        evo.frames_mean_p{c}(isnan(evo.frames_mean_p{c}))=0;
        evo.frames_mean_T{c}(isnan(evo.frames_mean_T{c}))=0;
        evo.cfc_samples_mean{c}(st) = mean(evo.cfc_samples{c}{st}(:));
        evo.cfc_samples_std{c}(st) = std(evo.cfc_samples{c}{st}(:));
        evo.cfc_samples_sem{c}(st) = std(evo.cfc_samples{c}{st}(:))/sqrt(length(evo.cfc_samples{c}{st}(:)));
        
    end
end

%end function
end





