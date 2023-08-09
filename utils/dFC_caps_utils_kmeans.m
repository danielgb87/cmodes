function dFC_caps_kmeans(data_chosen_parcel,nclust,dist,iter,reps,onlinePhase,start,opt)
% This function computes CAPs for each cluster size given.

%%%%%%%%%%%%%%%%%%%%%%%%

%now for every K, create a folder with the CAP mean maps, and T-maps.
Ncaps = nclust;
mkdir('CAP_results')
cd('CAP_results')
    for k = inputs.kmeans.Ncaps
        tmp = [];
        [tmp.frame_index, tmp.Centroids, tmp.sumd, tmp.Dist_to_centroid] = kmeans(data_chosen_seed1,k,...
            'Distance', inputs.kmeans.distance,...
            'MaxIter',inputs.kmeans.max_iter,...
            'OnlinePhase',inputs.kmeans.online_phase,...
            'Replicates', inputs.kmeans.Nreps,...
            'Start', inputs.kmeans.start,...
            'Options',inputs.kmeans.opts);

        % compute and map CAPs using data (not masked). The reduce results to
        % single format
        tmp = caps_compute_map_caps(data,tmp,k);

        tmp.frame_index = single(tmp.frame_index);
        tmp.cap_p_map_fdr = single(tmp.cap_p_map_fdr);

        %create directory for maps
        tmp_dir = (['CAPS_k' num2str(k)]);
        mkdir(tmp_dir)
        cd(tmp_dir)

        % Build maps
        for c = 1:k
            caps_create_nii(tmp.cap_mean_map(c,:), mask_info, ['CAP_' num2str(c) '_mean_BOLD']);
            caps_create_nii(tmp.cap_T_map(c,:), mask_info, ['CAP_' num2str(c) '_T_map']);
        end
        % recompute CAPs by adding garbage cluster
        tmp = caps_add_garbage_cluster(tmp,data);
        
        % re-map CAPs
        tmp_dir2 = (['CAPS_k+1_' num2str(k)]);
        mkdir(tmp_dir2)
        cd(tmp_dir2)
         for c = 1:k+1
            caps_create_nii(tmp.addk.cap_mean_map(c,:), mask_info, ['CAP_' num2str(c) '_mean_BOLD']);
            caps_create_nii(tmp.addk.cap_T_map(c,:), mask_info, ['CAP_' num2str(c) '_T_map']);
        end
        cd ..
        cd ..
        
        caps_results.(['CAPS_' num2str(k)]) = tmp;
        display(['done with k=' num2str(k)])
        
        
    end












end