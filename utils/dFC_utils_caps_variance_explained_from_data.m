function [r2,Ib,Iw,RD] = dFC_utils_caps_variance_explained_from_data(results, data)

[nobs, nvox]=size(data);
k=max(results.frame_index);
for c = 1:k
    
    % within cluster squared distances
    Nc(c) = nobs*results.occ_prob_mean(c);
    indc{c}=find(results.frame_index==c);
    
    for f = 1:length(indc{c})
        d2{c}(f) = (pdist2(results.Centroids(c,:),data(indc{c}(f),:),'correlation')).^2;
    end
    sumd2within(c) = sum(d2{c});
    
    Cw(c,:) = (Nc(c)/nobs).*results.Centroids(c,:);
end

% between cluster squared distances
Cbar=sum(Cw);         % weighted contributions to center of gravity.
for c = 1:k
    dist2Cbar(c) = Nc(c)*(pdist2(results.Centroids(c,:),Cbar,'correlation')).^2;
end


Iw = (1/nobs)*sum(sumd2within);   % within-cluster sum of
Ib = (1/nobs)*sum(dist2Cbar);   % cluster to Cg sum of distances

r2 = Ib/(Iw+Ib);

RD = Iw/Ib;

end %function