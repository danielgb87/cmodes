function caps_roi_kernel_at_CAP

% this function takes a CAPs analysis and a set of ROIs, and computes the
% time-locked kernel (mean +/- std) of the 15 previous and subsequent ROI
% BOLD values at points when the highest poeaks of a CAP occur.



%% 0. LOAD DATA, RESULTS, AND CFC TIMECOURSES. 

clc,clear
k = 6;
study_name = '16p11';
maindir = '/home/safaai/DanielG/caps_analysis/autism_study/16p11';
maindir_caps = [maindir '/runs_woCerebVent_centroids/k_analysis/cap_analysis_k' num2str(k)];
maindir_fc = [maindir '/seed_FC_at_CAP_analysis'];
roi_dir = '/home/safaai/DanielG/caps_analysis/autism_study/16p11/seed_FC_at_CAP_analysis/roi_kernels/';
cd(maindir)

% load data
load('mask_info_chd8_lowres')
load('16p11_data_wt_lowres'), data_1 = data; clear data
load('16p11_data_ht_lowres'), data_2 = data; clear data

% name datasets
dsname1 = 'WT';
dsname2 = 'HT';

% load CAP analysis results.
cd([maindir_caps '/caps_dynamics_k' num2str(k)])
load('caps_compare_analysis')
load('caps_compare_analysis_stats')
load('inputs_1')
load('inputs_2')

ord1 = caps_comp_analysis.cap_order_1;
ord2 = caps_comp_analysis.cap_order_2;

% optional, re-do the cpas_evo part

cap_evo_group_percent = 0.1; 
cap_evo_group_range = 30;
[caps_comp_analysis.evo_group.ds1] = caps_analysis_cap_evo_group(data_1, caps_comp_analysis.cap_CFC_analysis.ds1.CFC_peak_ind, caps_comp_analysis.cap_CFC_analysis.ds1.CFC_peak_amp, caps_comp_analysis.cap_occ.cap_occ_1_mean,cap_evo_group_range, caps_comp_analysis.cap_CFC_analysis.ds1.CFC, cap_evo_group_percent);
[caps_comp_analysis.evo_group.ds2] = caps_analysis_cap_evo_group(data_2, caps_comp_analysis.cap_CFC_analysis.ds2.CFC_peak_ind, caps_comp_analysis.cap_CFC_analysis.ds2.CFC_peak_amp, caps_comp_analysis.cap_occ.cap_occ_2_mean,cap_evo_group_range, caps_comp_analysis.cap_CFC_analysis.ds2.CFC, cap_evo_group_percent);



%% 1. Load ROIs, extract TS, and compute FC matrices
cd(roi_dir)
roi_files = dir('*.nii.gz');
nsubs1 = length(data_1);
nsubs2 = length(data_2);

for i = 1:length(roi_files)
    roi_names{i} = roi_files(i).name(1:end-7);
    roi_data_1{i} = caps_SB_extract_TS(data_1, [roi_files(i).name], roi_names{i}, mask_info.mask);
    roi_data_2{i} = caps_SB_extract_TS(data_2, [roi_files(i).name], roi_names{i}, mask_info.mask);
    %normalize signals
    roi_data_1{i}.TS = postproc_normalize(roi_data_1{i}.TS);
    roi_data_2{i}.TS = postproc_normalize(roi_data_2{i}.TS);
   
end
roi_list = roi_names;

%% 2. compute kernels at each CAP
% take the BOLD signal of the rois at the highest CAP occurrences (as
% included in caps_comp_analysis.evo_group.ds1.peaks_ind_final)

lags = -15:15;
for r = 1:length(roi_list)
    for c = 1:k
        for id = 1:length(caps_comp_analysis.evo_group.ds1.peaks_ind_final{ord1(c)})
            sub = caps_comp_analysis.evo_group.ds1.peaks_ind_final{ord1(c)}{id}(1);
            tt = caps_comp_analysis.evo_group.ds1.peaks_ind_final{ord1(c)}{id}(2);

            roi_kernel_at_cap_1{r,c}.roi_BOLD(id,:) = roi_data_1{r}.TS{sub}(tt+lags);
            
        end
        roi_kernel_at_cap_1{r,c}.BOLD_mean(:) = mean(roi_kernel_at_cap_1{r,c}.roi_BOLD,1);
        roi_kernel_at_cap_1{r,c}.BOLD_std(:) = std(roi_kernel_at_cap_1{r,c}.roi_BOLD,1);
        
        for id = 1:length(caps_comp_analysis.evo_group.ds2.peaks_ind_final{ord2(c)})
            sub = caps_comp_analysis.evo_group.ds2.peaks_ind_final{ord2(c)}{id}(1);
            tt = caps_comp_analysis.evo_group.ds2.peaks_ind_final{ord2(c)}{id}(2);
            
            roi_kernel_at_cap_2{r,c}.roi_BOLD(id,:) = roi_data_2{r}.TS{sub}(tt+lags);
            
        end
        roi_kernel_at_cap_2{r,c}.BOLD_mean(:) = mean(roi_kernel_at_cap_2{r,c}.roi_BOLD,1);
        roi_kernel_at_cap_2{r,c}.BOLD_std(:) = std(roi_kernel_at_cap_2{r,c}.roi_BOLD,1);
        
    end
end

% plot each activation function for each roi at each CAP
TR = 1.2;
for c = 1:k
    for r = 1:length(roi_names)
        
        figure
        shadedErrorBar(TR*lags,roi_kernel_at_cap_1{r,c}.BOLD_mean(:),roi_kernel_at_cap_1{r,c}.BOLD_std(:)/sqrt(size(roi_kernel_at_cap_1{r,c}.roi_BOLD,1)),'b',0.5)
        hold on
        shadedErrorBar(TR*lags,roi_kernel_at_cap_2{r,c}.BOLD_mean(:),roi_kernel_at_cap_2{r,c}.BOLD_std(:)/sqrt(size(roi_kernel_at_cap_2{r,c}.roi_BOLD,1)),'r',0.5)
        set(gcf,'Renderer','painters')
        set(gcf,'color',[0 0 0])
        %set(gcf, 'InvertHardCopy','off')
        ax=gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        box off
        set(gca,'color','none')
        hold on
        h1 = animatedline;
        h2 = animatedline;
        h1.Color = 'b';
        h2.Color = 'r';
        h1.Marker = 'o';
        h2.Marker = 'o';
        h1.MarkerFaceColor = 'b';
        h2.MarkerFaceColor = 'r';
        h1.MarkerSize = 12;
        h2.MarkerSize = 12;
        
        
        
        x = lags;
        y1 = double(roi_kernel_at_cap_1{r,c}.BOLD_mean(:));
        y2 = double(roi_kernel_at_cap_2{r,c}.BOLD_mean(:));
        
        a=tic; %start timer
        for t = 1:length(lags)
            addpoints(h1,x(t)*TR,y1(t))
            addpoints(h2,x(t)*TR,y2(t))
            b = toc(a); % check timer
            if b > 5 % it's the slowest mango can create videos
                drawnow
                a= tic;
            end
            %capture plot as img
            frame = getframe(gcf);
            img=frame2im(frame);
            [img,cmap]=rgb2ind(img,256);
            if t == 1
                imwrite(img,cmap,[roi_names{r} '_at_CAP' num2str(c) '.gif'], 'gif','LoopCount',Inf,'DelayTime',1);
            else
                imwrite(img,cmap,[roi_names{r} '_at_CAP' num2str(c) '.gif'], 'gif','WriteMode','append','DelayTime',1);
            end
            
        end
        
        close all
    end
end
    
    
    

end % function