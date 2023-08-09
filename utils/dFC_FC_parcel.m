function dFC_FC_parcel

% This is a function computes FC matrices from already extracted data (see dFC_data_prepare_parcel.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 7-03-2020
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% 1.a load data.
main = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/analysis_FC/FC_parcel_noScrub/';
main_data = '/media/DATA/dgutierrez/macaque_data/caps_concat_analysis_GM_downsampled/';
cd(main_data)
load('data_parcel_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
scrub = 0;
if scrub==1
    data_parcel = data_parcel_scrub;
    clear data_parcel_scrub
end
load('inputs_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
load('subject_list_mac_concat_2mm_woVent_WM_Cereb_BS.mat')
ndatasets = 3;
for nd = 1:ndatasets
    nsubs(nd) = length(data_parcel{nd});
end
nrois = inputs.Nrois;

%% 2. Compute FC for each dataset
cd(main)
for nd = 1:ndatasets
    for sub = 1:nsubs(nd)
        FC{nd}(sub,:,:) = corr(data_parcel{nd}{sub});
        FC_z{nd}(sub,:,:) = atanh(FC{nd}(sub,:,:));
        FC_z{nd}(isinf(FC_z{nd})|isnan(FC_z{nd})) = 0;
    end
    FC_mean{nd} = tanh(squeeze(mean(FC_z{nd},1)));
    for r1 = 1:nrois
        for r2 = 1:nrois
            [FC_h{nd}(r1,r2),FC_p{nd}(r1,r2),~,t] = ttest(spm_vec(FC_z{nd}(:,r1,r2)),0,'tail','both');
            FC_T{nd}(r1,r2) = t.tstat;
        end 
    end
    FC_T{nd}(isinf(FC_T{nd})|isnan(FC_T{nd})) = 0;
    FC_p{nd}(isinf(FC_p{nd})|isnan(FC_p{nd})) = 0;
    [~, fdr_crit_p(nd), ~]=fdr_bh_groppe(spm_vec(nonzeros(triu(FC_p{nd}))));

end

%% 3. compute FC differences between the datasets.
for nd1= 1:ndatasets
    for nd2=1:ndatasets
            for r1 = 1:nrois
                for r2 = 1:nrois
                    [FC_diff_h{nd1,nd2}(r1,r2), FC_diff_p{nd1,nd2}(r1,r2),~,t] = ttest2(spm_vec(FC_z{nd1}(:,r1,r2)),spm_vec(FC_z{nd2}(:,r1,r2)),'tail','both');
                    FC_diff_T{nd1,nd2}(r1,r2) = t.tstat;
                end 
            end
            FC_diff_T{nd1,nd2}(isinf(FC_diff_T{nd1,nd2})|isnan(FC_diff_T{nd1,nd2})) = 0;
            FC_diff_p{nd1,nd2}(isinf(FC_diff_p{nd1,nd2})|isnan(FC_diff_p{nd1,nd2})) = 0;
            [~, fdr_diff_crit_p(nd1,nd2), ~]=fdr_bh_groppe(spm_vec(nonzeros(triu(FC_diff_p{nd1,nd2}))));
    end
end

ncomp = nrois*(nrois-1)/2;
% FDR correct results
for nd1 = 1:ndatasets
    for nd2 = 1:ndatasets
        FC_diff_T_fdr{nd1,nd2}=FC_diff_T{nd1,nd2};
        for r1=1:nrois
            for r2 = 1:nrois
                if FC_diff_p{nd1,nd2}(r1,r2) > fdr_diff_crit_p(nd1,nd2)
                    FC_diff_T_fdr{nd1,nd2}(r1,r2)=0;
                end
            end
        end
    end
    FC_T_fdr{nd1}=FC_T{nd1};
    for r1=1:nrois
            for r2 = 1:nrois
                if FC_p{nd1}(r1,r2) > fdr_crit_p(nd1)
                    FC_T_fdr{nd1}(r1,r2)=0;
                end
            end
    end
end

% finally select the pairs of ROIs with the biggest between group
% differences by sorting the T-scores
kk=8;
for nd1 = 1:ndatasets 
    for nd2 = 1:ndatasets
        [diff_max20{nd1,nd2},diff_ind_max20{nd1,nd2}] = maxk(abs(FC_diff_T_fdr{nd1,nd2}),kk);
    end 
end

% look at the histogram of Tscores and select those differences above T=4
% that are reproducible
Tthr=2;
count=0;
for n = 1:nrois
    v1 = diff_max20{1,2}(:,n);
    v2 = diff_max20{1,3}(:,n);
    vi1 = diff_ind_max20{1,2}(:,n);
    vi2 = diff_ind_max20{1,3}(:,n);
        for m1 = 1:kk
            for m2 = 1:kk
                if vi1(m1)==vi2(m2) && v1(m1) > Tthr && v2(m2) > Tthr
                    count =count+1;
                    rp{count} = [n,vi1(m1)];
                    rp_Tmean(count) = mean([v1(m1),v2(m2)]);
                end
            end
        end
        rp_sum(n) = mean([sum(v1),sum(v2)]);
end
[rp_sort_T, rp_sort_ind] = sort(rp_Tmean, 'descend');               
[rpsum_sort_T, rpsum_sort_ind] = sort(rp_sum, 'descend');               

roi_pairs = rp(rp_sort_ind(1:2:end));
roi_pairs_Tmean = rp_sort_T(1:2:end);
for c = 1:length(roi_pairs_Tmean)
    if FC_diff_T_fdr{1,2}(roi_pairs{c}) < 0
        roi_pairs_Tmean(c) = -1*roi_pairs_Tmean(c);
    end
end 
%save into structure
FC_results.sort_FC_diff.diff_max20 = diff_max20;
FC_results.sort_FC_diff.diff_ind_max20 = diff_ind_max20;
FC_results.sort_FC_diff.roi_pairs = roi_pairs;
FC_results.sort_FC_diff.roi_pairs_Tmean = roi_pairs_Tmean;
FC_results.sort_FC_diff.roi_pairs_Tsum = rpsum_sort_T;
FC_results.sort_FC_diff.roi_pairs_Tsum_ind = rpsum_sort_ind;
            
%% 4. plot matrices
for nd1 = 1:ndatasets-1
    for nd2 = nd1+1:ndatasets
        figure('Position', [15,15,800,750])
        imagesc(FC_diff_T{nd1,nd2});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are      
        colorbar
        caxis([-5 5])
        fig.PaperPositionMode = 'auto';
        print(['FC_diff_T_ds' num2str(nd1) '-ds' num2str(nd2)], '-dpng','-r600')
        
        figure('Position', [15,15,800,750])
        imagesc(FC_diff_T_fdr{nd1,nd2});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are      
        colorbar
        caxis([-5 5])
        fig.PaperPositionMode = 'auto';
        print(['FC_diff_Tfdr_ds' num2str(nd1) '-ds' num2str(nd2)], '-dpng','-r600')
        close all
    end
end

for nd = 1:ndatasets
        figure('Position', [15,15,800,750])
        imagesc(FC_T{nd});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are      
        colorbar
        caxis([-8 8])
        fig.PaperPositionMode = 'auto';
        print(['FC_T_ds' num2str(nd)], '-dpng','-r600')
        
        figure('Position', [15,15,800,750])
        imagesc(FC_T_fdr{nd});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are      
        colorbar
        caxis([-8 8])
        fig.PaperPositionMode = 'auto';
        print(['FC_Tfdr_ds' num2str(nd)], '-dpng','-r600')
        
        figure('Position', [15,15,800,750])
        imagesc(FC_mean{nd});            % Create a colored plot of the matrix values
        colormap(othercolor('BuDRd_18'));  % Change the colormap to gray (so higher values are      
        colorbar
        caxis([-1 1])
        fig.PaperPositionMode = 'auto';
        print(['FC_mean_ds' num2str(nd)], '-dpng','-r600')
        close all 
end



%% 5. save results
FC_results.FC=FC;
FC_results.FC_diff_h=FC_diff_h;
FC_results.FC_diff_p=FC_diff_p;
FC_results.FC_diff_t=FC_diff_T;
FC_results.FC_diff_T_fdr=FC_diff_T_fdr;
FC_results.FC_h=FC_h;
FC_results.FC_mean = FC_mean;
FC_results.FC_p=FC_p;
FC_results.FC_T=FC_T;
FC_results.FC_T_fdr=FC_T_fdr;
FC_results.FC_z=FC_z;
FC_results.fdr_crit_p=fdr_crit_p;
FC_results.fdr_diff_crit_p = fdr_diff_crit_p;
FC_results.ncomp=ncomp;
FC_results.ndatasets=ndatasets;
FC_results.nrois = nrois;
save(['FC_results'] ,'FC_results','-v7.3')



% end main function.
end



