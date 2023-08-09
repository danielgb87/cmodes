function dFC_utils_caps_RSN_in_caps

% this functiona takes a set of resting-state networkswithin the same space
% of CAPs results, and models each CAP mean map as a linear combination of
% these binary RSN maps.


clc, clear

%% 1. organize and load data
% define a directory to save, and load the CAPs results and rsn maps with
% id
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_ds/cap_maps_2mm_k8/rsn_in_caps';
rsn_dir = '/media/DATA/dgutierrez/macaque_data/Ting_parcellation_2019/rsns';
load ('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/mask_info_mac_concat_2mm_woVent_WM_Cereb_BS')
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_ds/cap_maps_2mm_k8/cap_results_k8')
k=8;
nds = 3;
cd(rsn_dir)
rsn_file = dir('*.nii*');
for m=1:length(rsn_file)   
    % Go to the subjects preprocessing dir 
        temp = spm_vol(rsn_file(m).name); 
        [~, rsn_name{m}, ~] = fileparts(rsn_file(m).name);
        data1 = spm_vec(spm_read_vols(temp));      
        data = data1(spm_vec(mask_info.mask)==1,:);
        rsn(m,:) = data; 
        clear data1 temp
end
nr = length(rsn_file);


%% 2. for each CAP, make a GLM
cd(main)
for nd = 1:nds
    caps{nd} = results_map.map_vw_ds{nd}.cap_T_map;
    caps_p{nd} = results_map.map_vw_ds{nd}.cap_p_map;
    for c = 1:k
%         [glm_cap{nd,c},dev{nd,c},stats{nd,c}] = glmfit(rsn',spm_vec(caps{nd}(c,:)));
        cmap = spm_vec(caps{nd}(c,:)); cmap_thr=cmap;
        pmap = spm_vec(caps_p{nd}(c,:));
        cmap_thr(pmap > 0.01) = 0;
        mdl{nd,c} = fitglm(rsn',cmap,'Intercept',false);
        mdl_thr{nd,c} = fitglm(rsn',cmap_thr,'Intercept',false);
        
        r2(nd,c) = mdl{nd,c}.Rsquared.Adjusted;
        r2_thr(nd,c) = mdl_thr{nd,c}.Rsquared.Adjusted;
        for m = 1:nr
            rsn_coef{nd}(c,m) = mdl{nd,c}.Coefficients{['x' num2str(m)],'Estimate'};
            rsn_p{nd}(c,m) = mdl{nd,c}.Coefficients{['x' num2str(m)],'pValue'};
            rsn_T{nd}(c,m) = mdl{nd,c}.Coefficients{['x' num2str(m)],'tStat'};
            
            rsn_coef_thr{nd}(c,m) = mdl_thr{nd,c}.Coefficients{['x' num2str(m)],'Estimate'};
            rsn_p_thr{nd}(c,m) = mdl_thr{nd,c}.Coefficients{['x' num2str(m)],'pValue'};
            rsn_T_thr{nd}(c,m) = mdl_thr{nd,c}.Coefficients{['x' num2str(m)],'tStat'};
        end
    end
    
end

% create significance matrices
for nd = 1:nds    
    rsn_h_thr{nd} = zeros(k,nr);
    rsn_h_thr{nd}(rsn_p_thr{nd} < 0.01)=1;
    rsn_h{nd} = zeros(k,nr);
    rsn_h{nd}(rsn_p_thr{nd} < 0.01)=1;
end

%% 3. create the RSN to cap matrices and plot
for nd = 1:nds
        mat = rsn_coef_thr{nd}.*rsn_h_thr{nd};
        figure
        imagesc(mat);    
        colormap(othercolor('BuDRd_18'))
        hold on;
        axis image
        axis ij
        colorbar
        caxis([-10 10]) 
        pbaspect([k nr 1])
        set(gca,'xtick', linspace(0,k+0.5,k+0.5), 'ytick', linspace(0,nr+0.5,nr+0.5));
        set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
        for c = 1:size(mat,1)
            cnames{c} = ['CAP ' num2str(c)];
        end
        for m = 1:size(mat,2)
            rnames{m} = rsn_name{m}(1:end-4);
        end
        set(gca, 'XTick', 1:size(mat,2), ...                             % Change the axes tick marks
                 'XTickLabel', rnames, ...  %   and tick labels
                 'YTick', 1:size(mat,1), ...
                 'YTickLabel', cnames, ...
                 'TickLength', [0 0]);
                  xtickangle(45)
        fig.PaperPositionMode = 'auto';
        print(['RSNs_in_caps_ds' num2str(nd)], '-dpng','-r600')
        close all

end
    








%end function
end