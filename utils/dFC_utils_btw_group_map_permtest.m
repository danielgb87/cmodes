function [rcrit, rvec, pval]=dFC_utils_btw_group_map_permtest(data1, data2, nsel1, nsel2, rtrue, nperm)

% This script takes two data matrices (ntimepoints x ndimensions, and
% randomly selects "nsel1" and "nsel2" timepoints and compares the
% similarity of the surrogate means (in time) of each data matrix.
for p = 1:nperm
        c1=randperm(size(data1,1),nsel1);
        c2=randperm(size(data2,1),nsel2);
        dsurr1=data1(c1,:);
        dsurr2=data2(c2,:);
        map1=mean(dsurr1,1);
        map2=mean(dsurr2,1);
        rvec(p)=corr(map1',map2');
end
rcrit=max(rvec);
pval=sum(rvec>rtrue)/nperm;
    
    
    








end