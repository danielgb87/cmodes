function [T] = dFC_utils_FC_group_ttest(nvar, FC)

%stat tests
FC_group_T=zeros(nvar,nvar);
FC_group_h=zeros(nvar,nvar);
FC_group_p=zeros(nvar,nvar);
for v1 = 1:nvar-1
    parfor v2 = v1+1:nvar
        [FC_group_h(v1,v2),FC_group_p(v1,v2),~,t] = ttest(spm_vec(FC(:,v1,v2)),0,'tail','both');
        FC_group_T(v1,v2) = t.tstat;
    end    
end
%copy diagonal matrix and eliminate IsInf and IsNaN
FC_group_p(isnan(FC_group_p)) = 1;
FC_group_h(isnan(FC_group_p)) = 0;
FC_group_T(isinf(FC_group_T)|isnan(FC_group_T)) = 0;
T.FC_group_T=FC_group_T + FC_group_T';
T.FC_group_p=FC_group_p + FC_group_p';
T.FC_group_h=FC_group_h + FC_group_h';

%FDR correct
[~, FC_group_p_crit, ~]=fdr_bh_groppe(spm_vec(nonzeros(triu(FC_group_p))));
FC_group_T_fdr=FC_group_T;
for v1=1:nvar-1
    for v2 = v1+1:nvar
        if FC_group_p(v1,v2) > FC_group_p_crit
            FC_group_T_fdr(v1,v2)=0;
        end
    end
end
T.FC_group_T_fdr=FC_group_T_fdr+FC_group_T_fdr';
T.FC_group_p_crit=FC_group_p_crit;


end %function