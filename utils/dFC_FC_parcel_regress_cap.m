function dFC_FC_parcel_regress_cap

% This script takes the FC info and CAP info, and recomputes FC after
% regressing each cap. Make sure to order the CAPs before starting in the editing section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 7-03-2020
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% 1.a load data.
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC/FC_parcel/';
main_data = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/';
cd(main_data)
load('data_parcel_scrub_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
scrub = 1;
if scrub==1
    data_parcel = data_parcel_scrub;
    clear data_parcel_scrub
end
analisys_id = 'FC_regCap_parcel';
load('/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_parcel/cap_maps_2mm_k8/caps_dynamics_k8/caps_dyn')
load('inputs_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('subject_list_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
cap_ord = [1 3 2 6 4 7 5 8];
k=8;
ndatasets = 3;
for nd = 1:ndatasets
    nsubs(nd) = length(data_parcel{nd});
end
nrois = inputs.Nrois;

%% 2. compute FC after regressing each CAP
cd(main)
for c = 1:k
    for nd = 1:ndatasets
        cfc = caps_dyn.cfc.cfc{nd};
        for sub = 1:nsubs(nd)
            for n = 1:nrois
                    [~,~,data_reg{nd,c}{sub}(:,n)] = regress(spm_vec(data_parcel{nd}{sub}(:,n)),[spm_vec(ones(length(cfc.cfc_norm{sub,c}),1)) spm_vec(cfc.cfc_norm{sub,c})]);
            end
        end
    end
end
for nd = 1:ndatasets
        for sub = 1:nsubs(nd)
            for n = 1:nrois
                    [~,~,data_reg{nd,k+1}{sub}(:,n)] = regress(spm_vec(data_parcel{nd}{sub}(:,n)),[spm_vec(ones(length(caps_dyn.gs.gs_norm{nd}{sub}),1)) spm_vec(caps_dyn.gs.gs_norm{nd}{sub})]);
            end
        end
end


for c = 1:k+1
    for nd = 1:ndatasets
        for sub = 1:nsubs(nd)
            
            FC_reg{nd,c}(sub,:,:) = corr(data_reg{nd,c}{sub});
            FC_z_reg{nd,c}(sub,:,:) = atanh(FC_reg{nd,c}(sub,:,:));
            FC_z_reg{nd,c}(isinf(FC_z_reg{nd,c})|isnan(FC_z_reg{nd,c})) = 0;       
        end
        FC_mean_reg{nd,c} = tanh(squeeze(mean(FC_z_reg{nd,c},1)));
        for r1 = 1:nrois
            for r2 = 1:nrois
                [FC_reg_h{nd,c}(r1,r2),FC_p_reg{nd,c}(r1,r2),~,t] = ttest(spm_vec(FC_z_reg{nd,c}(:,r1,r2)),0,'tail','both');
                FC_T_reg{nd,c}(r1,r2) = t.tstat;
            end 
        end
        FC_T_reg{nd,c}(isinf(FC_T_reg{nd,c})|isnan(FC_T_reg{nd,c})) = 0;
        FC_p_reg{nd,c}(isinf(FC_p_reg{nd,c})|isnan(FC_p_reg{nd,c})) = 0;
        [~, fdr_crit_p(nd,c), ~]=fdr_bh_groppe(spm_vec(nonzeros(triu(FC_p_reg{nd,c}))));
        FC_T_fdr_reg{nd,c}=FC_T_reg{nd};
        for r1=1:nrois
                for r2 = 1:nrois
                    if FC_p_reg{nd,c}(r1,r2) > fdr_crit_p(nd,c)
                        FC_T_fdr_reg{nd,c}(r1,r2)=0;
                    end
                end
        end
    end
    
    % compute differences
    
    for nd1= 1:ndatasets
        for nd2=1:ndatasets
            for r1 = 1:nrois
                for r2 = 1:nrois
                    [FC_diff_h_reg{c}{nd1,nd2}(r1,r2), FC_diff_p_reg{c}{nd1,nd2}(r1,r2),~,t] = ttest2(spm_vec(FC_z_reg{nd1,c}(:,r1,r2)),spm_vec(FC_z_reg{nd2,c}(:,r1,r2)),'tail','both');
                    FC_diff_T_reg{c}{nd1,nd2}(r1,r2) = t.tstat;
                end
            end
            FC_diff_T_reg{c}{nd1,nd2}(isinf(FC_diff_T_reg{c}{nd1,nd2})|isnan(FC_diff_T_reg{c}{nd1,nd2})) = 0;
            FC_diff_p_reg{c}{nd1,nd2}(isinf(FC_diff_p_reg{c}{nd1,nd2})|isnan(FC_diff_p_reg{c}{nd1,nd2})) = 0;
            [~, fdr_diff_crit_p_reg{c}(nd1,nd2), ~]=fdr_bh_groppe(spm_vec(nonzeros(triu(FC_diff_p_reg{c}{nd1,nd2}))));
            FC_diff_T_fdr_reg{c}{nd1,nd2}=FC_diff_T_reg{c}{nd1,nd2};
            for r1=1:nrois
                for r2 = 1:nrois
                    if FC_diff_p_reg{c}{nd1,nd2}(r1,r2) > fdr_diff_crit_p_reg{c}(nd1,nd2)
                        FC_diff_T_fdr_reg{c}{nd1,nd2}(r1,r2)=0;
                    end
                end
            end
        end
    end
end


%% 3. plot results
for c = 1:k+1
    for nd1 = 1:ndatasets-1
        for nd2 = nd1+1:ndatasets
            figure('Position', [15,15,800,750])
            imagesc(FC_diff_T_reg{c}{nd1,nd2});            % Create a colored plot of the matrix values
            colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are
            colorbar
            caxis([-5 5])
            fig.PaperPositionMode = 'auto';
            print(['FC_diff_T_reg_C' num2str(c) '_ds' num2str(nd1) '-ds' num2str(nd2)], '-dpng','-r600')
            
            figure('Position', [15,15,800,750])
            imagesc(FC_diff_T_fdr_reg{c}{nd1,nd2});            % Create a colored plot of the matrix values
            colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are
            colorbar
            caxis([-5 5])
            fig.PaperPositionMode = 'auto';
            print(['FC_diff_Tfdr_reg_C' num2str(c) '_ds' num2str(nd2)], '-dpng','-r600')
            close all
        end
    end
    
    for nd = 1:ndatasets
        figure('Position', [15,15,800,750])
        imagesc(FC_T_reg{nd,c});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are
        colorbar
        caxis([-8 8])
        fig.PaperPositionMode = 'auto';
        print(['FC_T_reg_C' num2str(c) '_ds' num2str(nd)], '-dpng','-r600')
        
        figure('Position', [15,15,800,750])
        imagesc(FC_T_fdr_reg{nd,c});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are
        colorbar
        caxis([-8 8])
        fig.PaperPositionMode = 'auto';
        print(['FC_T_reg_C' num2str(c) '_ds' num2str(nd)], '-dpng','-r600')
        
        figure('Position', [15,15,800,750])
        imagesc(FC_mean_reg{nd,c});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are
        colorbar
        caxis([-1 1])
        fig.PaperPositionMode = 'auto';
        print(['FC_mean_reg_C' num2str(c) '_ds' num2str(nd)], '-dpng','-r600')
        close all
    end
    
end


    
    


end % function