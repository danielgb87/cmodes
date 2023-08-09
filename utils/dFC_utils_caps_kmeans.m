function [caps_results] = dFC_utils_caps_kmeans(data,nsubs,nobs,nclust,dist,iter,reps,onlinePhase,start,opt)
% This function computes CAPs for each cluster size given.

%%%%%%%%%%%%%%%%%%%%%%%%

%now for every K, create a folder with the CAP mean maps, and T-maps.
Ncaps = nclust;

for k = Ncaps
    tmp = [];
    [tmp.frame_index, tmp.Centroids, tmp.sumd, tmp.Dist_to_centroid] = kmeans(data,k,...
        'Distance', dist,...
        'MaxIter',iter,...
        'OnlinePhase',onlinePhase,...
        'Replicates', reps,...
        'Start', start,...
        'Options',opt);
    
    tmp.frame_index = single(tmp.frame_index);
        
    %compute between CAP similarity:
    tmp.btw_cap_sim = corr(tmp.Centroids');

    
    % put results in single-suject form
        % count CAP occurrences for each subject, and compute their proportion.
        t=1;
        for sub = 1:nsubs
            tmp.frame_ind_sub{sub}(:,1) = tmp.frame_index(t:t+nobs(sub)-1);
            t = t+nobs(sub);

            for c=1:k
                tmp.occ_prob_sub(sub,c) = sum(tmp.frame_ind_sub{sub}==c)/nobs(sub);
                tmp.occ_sub(sub,c) = sum(tmp.frame_ind_sub{sub}==c);

            end

        end

        % compute mean and std.
        for c=1:k
            tmp.occ_prob_mean(c) = mean(tmp.occ_prob_sub(:,c));
            tmp.occ_prob_std(c) = std(tmp.occ_prob_sub(:,c));

            tmp.occ_mean(c) = mean(tmp.occ_sub(:,c));
            tmp.occ_std(c) = std(tmp.occ_sub(:,c));
        end
        
    caps_results.(['caps_' num2str(k)]) = tmp;
    display(['done with k=' num2str(k)])
    
    
end

end

