function    [seed_data] = dFC_utils_SBA_ts_extract(data, seed_name, mask)

seed_full_img = spm_read_vols(spm_vol(seed_name));
% crop seed
vec = spm_vec(seed_full_img);
vec = vec(spm_vec(mask)>=1);
seed_ind = find(vec>=1);
for sub = 1:length(data)
    TS{sub,1} = mean(data{sub}(:,seed_ind),2);
end
seed_data.img = seed_full_img;
seed_data.name = seed_name;
seed_data.TS = TS;
% end funtion
end







