function tmp = dFC_utils_caps_mapping_ss(data,tmp,k)

for sub = 1:length(data)
    ind = tmp.frame_ind_sub{sub};
    for c = 1:k
        tmp.cap_mean_map_ss{sub}(c,:)  = mean(data{sub}(ind==c,:),1);
        
        % Perform t-test voxel-wise.
        data_temp = data{sub}(ind==c,:);
        parfor v = 1:size(data{sub},2)
            [~, pmap(v,1), ~, stats] = ttest(data_temp(:,v));
            Tmap(c,v) = stats.tstat;
        end
        
        
    end
    tmp.cap_p_map_ss{sub} = pmap;
    tmp.cap_T_map_ss{sub} = Tmap;
    tmp.cap_p_map_ss{sub}(isnan(tmp.cap_p_map_ss{sub}))=0;
    tmp.cap_T_map_ss{sub}(isnan(tmp.cap_T_map_ss{sub}))=0;
end


end