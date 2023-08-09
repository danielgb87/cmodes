function [Params] = dFC_utils_SBA_FCparams(TS,data,CMap)

% This function computes correlation map to volume averaging similarity
% using different thresholds. It was modified from the 'Comp_Params.m'
% function of the CAPsTOOLBOX(C) [1]
%
% [1]: Amico, Enrico, et al. "Posterior Cingulate Cortex-Related
% Co-Activation Patterns: A Resting State fMRI Study in Propofol-Induced
% Loss of Consciousness." PloS one 9.6 (2014): e100012.
%
% Written by EA Sep 2014.
%__________________________________________________________________________
% License

%__________________________________________________________________________

Params ={};

% For each subject, using the dataset (Nvols x Nvox), the seed time-course,
% and the correlation map, compute by steps the correlation to
% volume-averaging similarity.
for sub=1:length(data)
    CorrMap = CMap{sub};
    dat     = data{sub};
    dat(isnan(dat))=0;
    
    ts = TS{sub};                      % call seed Time-series
    Nvols = max(size(ts));                  % number of time-points
    
    % take maximum and minimum values of the seed-timecourse and a
    % range of value ths to test.
    [maxTS indmax] = max(ts);
    [minTS indmin] = min(ts);
    Start = maxTS;
    Step = (maxTS-minTS)/50;
    count = 1;
    
    for thr = Start:-Step:minTS                                 % for each thr value
        index1 = find(ts >= thr);                      % create a vector of indexes (TS >=thr)
        IndStruct{sub,count} = index1;                        % save the ind vector)
        ActiMap1 =zeros(1,size(dat,2));                % create activation map
        tmp = mean(dat(index1,:),1);                   % take the active time-points' frames and average them
        ActiMap1 = tmp;   % feed the mean activation map.
        rate1(sub,count) = (1-(length(index1)/Nvols))*100;     % compute the proportion of active time-points (rate)
        SpatCorr1(sub,count) = corr(CorrMap',ActiMap1');      %Spatial corr between mean activation and correlation map. needs transpose for this kind of arrays
        count = count+1;
    end
end

Params.frame_indexes = IndStruct;
Params.percentile = rate1;
Params.SpatCorr = SpatCorr1;


return;
