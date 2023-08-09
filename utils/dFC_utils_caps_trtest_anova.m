function dFC_caps_utils_trtest_anova_tukey


clc, clear

main = '/media/DATA/dgutierrez/human_data/caps_hnu_10sessions_GM_concat/analysis_parcel';
cd(main)
load('caps_props_concat_parcel.mat')
ncaps = 2:20;
nds = 10;
nsubs = 30;

for k = 1:length(ncaps)
    for n = 1:nds
        for sub = 1:nsubs
            for c = 1:ncaps(k)
                occ{k}{c}(sub,n) = caps_props.cap_occ_group{k,n}(sub,c);
            end
        end
        gname{n} = ['s' num2str(n)];
    end
end

for k = 1:length(ncaps)
    for c = 1:ncaps(k)
        [anv_p{k}(c)]= anova1(occ{k}{c},gname,'off');
        if anv_p{k}(c) < 0.05/ncaps(k)
            anv_p_flag_bonf(k) = 1;
        else
            anv_p_flag_bonf(k) = 0;
        end
        
    end
    [anv_fdr_h{k}(:) anv_p_crit(k) adj_p{k}(:)]=fdr_bh_groppe(anv_p{k},0.05);
    anv_fdr_hsum(k) = sum(anv_fdr_h{k}(:));
    for c = 1:ncaps(k)
        if anv_p{k}(c) < anv_p_crit(k)
            anv_p_flag_fdr(k) = 1;
        else
            anv_p_flag_fdr(k) = 0;
        end
    end
    
end

end %function






end % function