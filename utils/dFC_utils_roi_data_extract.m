function [roi_data] = dFC_utils_roi_data_extract(data, roi_path, roi_name, mask)

% This function takes the data of a subject in matrix form (within a mask)
% and extracts the signal from a given ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Daniel Gutierrez-Barragan 2020, V1 25-02-2020

roi_full_img = spm_read_vols(spm_vol(roi_path));
% crop roi
vec = spm_vec(roi_full_img);
vec = vec(spm_vec(mask)>=1);
roi_ind = find(vec>=1);
for sub = 1:length(data)
    TS{sub,1} = mean(data{sub}(:,roi_ind),2);
end

roi_data.img = roi_full_img;
roi_data.name = roi_name;
roi_data.TS = TS;


% end funtion
end







