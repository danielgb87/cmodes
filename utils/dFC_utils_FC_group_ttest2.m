function [T] = dFC_utils_FC_group_ttest2(nvar, FC1, FC2)

%stat tests
FC_group_diff_T=zeros(nvar,nvar);
FC_group_diff_h=zeros(nvar,nvar);
FC_group_diff_p=zeros(nvar,nvar);
for v1 = 1:nvar-1
    parfor v2 = v1+1:nvar
        [FC_group_diff_h(v1,v2),FC_group_diff_p(v1,v2),~,t] = ttest2(spm_vec(FC1(:,v1,v2)),spm_vec(FC2(:,v1,v2)),'tail','both');
        FC_group_diff_T(v1,v2) = t.tstat;
    end    
end
%copy diagonal matrix and eliminate IsInf and IsNaN
FC_group_diff_p(isnan(FC_group_diff_p)) = 1;
FC_group_diff_h(isnan(FC_group_diff_p)) = 0;
FC_group_diff_T(isinf(FC_group_diff_T)|isnan(FC_group_diff_T)) = 0;
T.FC_group_diff_T=FC_group_diff_T + FC_group_diff_T';
T.FC_group_diff_p=FC_group_diff_p + FC_group_diff_p';
T.FC_group_diff_h=FC_group_diff_h + FC_group_diff_h';

%FDR correct
[~, FC_group_diff_p_crit, ~]=fdr_bh_groppe(spm_vec(nonzeros(triu(FC_group_diff_p))));
FC_group_diff_T_fdr=FC_group_diff_T;
for v1=1:nvar-1
    for v2 = v1+1:nvar
        if FC_group_diff_p(v1,v2) > FC_group_diff_p_crit
            FC_group_diff_T_fdr(v1,v2)=0;
        end
    end
end
T.FC_group_diff_T_fdr=FC_group_diff_T_fdr + FC_group_diff_T_fdr';
T.FC_group_diff_p_crit=FC_group_diff_p_crit;


end %function