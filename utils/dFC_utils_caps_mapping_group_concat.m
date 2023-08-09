function tmp = dFC_utils_caps_mapping_group_concat(data,tmp,k)

ind = tmp.frame_index;
for c = 1:k
    tmp.cap_mean_map(c,:)  = mean(data(ind==c,:),1);
    
    % Perform t-test voxel-wise.
    data_temp = data(ind==c,:);
    parfor v = 1:size(data,2)
        [~, cap_p_map(v,1), ~, stats] = ttest(data_temp(:,v));
        cap_T_map(v,1) = stats.tstat;
    end
    cap_p_map(isnan(cap_p_map))=0;
    cap_T_map(isnan(cap_T_map))=0; 
    tmp.cap_p_map(c,:) = cap_p_map;
    tmp.cap_T_map(c,:) = cap_T_map;
      
    tmp.occ_prob_group(c)  = sum(ind == c)/length(ind);
    tmp.cap_consistency(c) = mean(corr(tmp.cap_mean_map(c,:)',data(ind==c,:)'));
        
end



end