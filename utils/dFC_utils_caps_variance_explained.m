function [Iw,Ib,VE,RD] = dFC_utils_caps_variance_explained(results)

nobs = length(results.frame_index);
[k, nvox] = size(results.Centroids);

wtn_clust_dist = zeros(size(results.Dist_to_centroid));
for c = 1:k
    for v = 1:nobs
        if results.frame_index(v) == c
            wtn_clust_dist2(v,c)=results.Dist_to_centroid(v,c)^2;
        end
    end
    wtn_clust_sumdist2(c) = sum(wtn_clust_dist2(:,c));
end
        
        
Iw = (1/nobs)*sum(wtn_clust_sumdist2);

for c = 1:k
    Nc(c) = sum(results.frame_index==c);
    Cbar(c,:) = (Nc(c)/nobs)*results.Centroids(c,:);
end

Cbar_w = sum(Cbar);

for c = 1:k
    d2Cbar(c) = Nc(c)*(pdist2(results.Centroids(c,:),Cbar_w,'correlation')).^2;
end

Ib = (1/nobs)*sum(d2Cbar);

RD = Iw/Ib;
VE = Ib/(Ib+Iw);


end %function

    