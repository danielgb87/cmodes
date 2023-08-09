function     [selected_frames, frame_list] = dFC_utils_SBA_seed_select_frames(TS, thr,params)

Perc = params.percentile;    % extract the threshold rates of each subject
Index = params.frame_indexes;    % extract the indices of supra-thr seed - activity for each subject and thr level

t = thr;
for sub=1:length(TS)
    selected_frames{sub} = zeros(length(TS{sub,1}),1);
    
    index = find(Perc(sub,:)>=(100-t));      % select the time indices of supra-thr seed-activity.
    IndFrames = Index{sub,index(end)};         % localize the frames
    selected_frames{sub}(IndFrames) = 1;
    frame_list{sub}(:) = find(selected_frames{sub}==1);
    
end


end

