function  [corrMap, corrMap_mean, corrMapT] = dFC_utils_SBA_FCcomp(data,TS)


for sub = 1:length(data)    
    corrMap{sub,1} = corr(TS{sub,1}, data{sub,1});
    corrMap{sub,1}(isnan(corrMap{sub,1})) = 0;    
end

% compute and map t-map
data_temp = cell2mat(corrMap);
corrMap_mean = mean(data_temp,1);

parfor v = 1:size(data_temp,2)
    [~, corr_map_p(1,v), ~, stats] = ttest(data_temp(:,v));
    corrMapT(1,v) = stats.tstat;
end

corr_map_p(isnan(corr_map_p))=0;
corrMapT(isnan(corrMapT))=0;


return;