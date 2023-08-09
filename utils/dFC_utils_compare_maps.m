function dFC_utils_compare_maps

% this script assign maps from one group, extracting form their .nii files
% to another group. Be sure to have the group with the biggest amount of
% element as the second group.
clc,clear
main = '/media/DATA/dgutierrez/human_data/caps_human/msc_5sessions_craddock950_ppTingDGB/analysis_parcel_scrub/cap_maps_3mm_k8/';
comp_name = 'MSCvsMRC_MSC_MSCMask';
dirmaps1 = '/media/DATA/dgutierrez/human_data/caps_human/msc_5sessions_craddock950_ppTingDGB/analysis_parcel_scrub/cap_maps_3mm_k8/';
dirmaps2 = '/media/DATA/dgutierrez/human_data_mrc_fabio/TD-males/K8/';
%define the number of maps in each group
nmaps1 = 8; 
nmaps2 = 8;
%load a common mask to compare only within mask voxels
mask_file = '/media/DATA/dgutierrez/human_data/masks_templates_parcels/fsl_masks/Mask_Craddock_parcel950_rpi_3mm.nii.gz';
mask = single(spm_read_vols(spm_vol(mask_file)));
% Check indexes with non-zero values
[mask_ind(1,:), mask_ind(2,:), mask_ind(3,:)] = ind2sub(size(mask),find(mask ==1));
mask_info.file = mask_file;
mask_info.mask = mask;
mask_info.mask_index = mask_ind;

% load caps from both groups
cd(dirmaps1)
cmap1 = cell(nmaps1,1);
mlist1 = dir('*concat_mean_cap*');
for c = 1:nmaps1
    temp = spm_vol(mlist1(c).name);
    data1 = spm_vec(spm_read_vols(temp));
    cmap1{c} = data1(spm_vec(mask)==1,:);
    cmap1{c,1} = cmap1{c,1}';
    clear data1 temp
end

cd(dirmaps2)
cmap2 = cell(nmaps2,1);
mlist2 = dir('*cap_mean_All_Sbj*');
for c = 1:nmaps2
    temp = spm_vol(mlist2(c).name);
    data2 = spm_vec(spm_read_vols(temp));
    
    cmap2{c} = data2(spm_vec(mask)==1,:);
    cmap2{c,1} = cmap2{c,1}';
    clear data2 temp
end

%% compare mean maps (match then compare) with the hungarian algorithm
% the munkres_HA function is needed.
for i = 1:nmaps1
    for j = 1:nmaps2
        Dmat(i,j) = pdist([(spm_vec(cmap1{i})'); (spm_vec(cmap2{j})')]);
    end
end
[x_ind, cost] = munkres_HA(Dmat);
%if you want to change order
% x_ind = [7 5 6 1];
%resolve the missing maps
nmiss=setdiff(1:nmaps2,x_ind);
x_ind_full = [x_ind, nmiss];

% Compute pairwise similarity between matched CAPs from each ds.
for c = 1:nmaps1
    btw_ds_sim(c) = corr(spm_vec(cmap1{c}),spm_vec(cmap2{x_ind(c)}));
    for c2 = 1:nmaps2
        corr_mat(c,c2) = corr(spm_vec(cmap1{c}),spm_vec(cmap2{x_ind_full(c2)}));
    end
end

%plot results
cd(main)
for c = 1:nmaps1
    cnames1{c} = ['map ' num2str(c)];
end
for c = 1:nmaps2
    cnames2{c} = ['map ' num2str(x_ind_full(c))];
end
figure
imagesc(corr_mat);            % Create a colored plot of the matrix values
colormap(jet);  % Change the colormap to gray (so higher values are

set(gca, 'XTick', 1:size(corr_mat,2), ...                             % Change the axes tick marks
         'XTickLabel', cnames2, ...  %   and tick labels
         'YTick', 1:size(corr_mat,1), ...
         'YTickLabel', cnames1, ...
         'TickLength', [0 0]);
xtickangle(45)
colorbar
caxis([-1 1])  
print(['btw_ds_corr_' comp_name], '-dpng','-r600')

close all


end %function