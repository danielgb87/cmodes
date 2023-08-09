function caps_parcel_cluster_centroid_stability_evolution(Max_caps)

% this function take results from clusterings with diverse k, and compares
% the stability of the clusters as k rises by comparing them to the
% previous k (k-1)


%% 1. load the results
load('caps_results_All_Sbj.mat')
ncaps = 2:Max_caps;


%% 2. for each k, assign the best fitting maps of k, to theprevious (k-1) caps
cap_maps{1} = caps_results{1}.cap_mean_map;
cap_corr = cell(length(ncaps),1);
cap_corr{1} =[0,0];
cap_ord{1} = [1,2];

for k =2:length(ncaps)   
    tmp1 = caps_results{k-1}.cap_mean_map;  
    tmp2 = caps_results{k}.cap_mean_map;
    % assign the caps from k to the k-1 clustering
    for i = 1:ncaps(k)-1
        for j = 1:ncaps(k)
            Dmat_m{k}(i,j) = pdist([spm_vec(tmp1(cap_ord{k-1}(i),:))'; spm_vec(tmp2(j,:))']);
        end
    end
    [cap_ord{k}, cost_m{k}] = munkres_HA(Dmat_m{k});
    % identify the "new CAP" in k
    v = 1:ncaps(k); 
    for c = 1:ncaps(k)
       if ismember(cap_ord{k},v(c)) ==0
          cap_ord{k}(k+1) = c;
       end
    end
    
    %compute the similarities between matched caps
    for c = 1:ncaps(k-1)
        cap_corr{k}(c) = corr(spm_vec(tmp1(cap_ord{k-1}(c),:)), spm_vec(tmp2(cap_ord{k}(c),:)));
    end
    cap_corr{k}(end+1) = 0;
     
end

% organize evolution vectors
cap_evo = cell(length(ncaps),1);

for k = 1:length(ncaps)
    for kk = k:length(ncaps)
        cap_evo{k}(kk) = cap_corr{kk}(k);
    end
end

%% 3. plot results
figure
fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 35;
fig.Position(4) = 10;
fig.Position(1) = 0;
fig.Position(2) = 0;
plot(2:length(ncaps)+1,cap_evo{1}(1:length(ncaps)));
hold on
for k = 2:length(cap_evo)
    plot(k:length(ncaps)+1,cap_evo{k}(k-1:length(ncaps)))
    hold on
end
hold on
xlabel('K','fontsize',12)
ylabel('CAP stability','fontsize',12)
xlim([0 ncaps(end)+1])
ylim([0 1])
title('Evolution of CAP stability')


fig.PaperPositionMode = 'auto';
print(['cap_evolution'], '-dpng','-r600');
    

end %function