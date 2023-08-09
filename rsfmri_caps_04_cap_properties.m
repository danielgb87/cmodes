function caps_04_dFC_caps_ds_properties


% this function performs various analyses to compare datasets. The results chosen are the ones from the run with the highest variance explained. It is
% intended for data that was concatenated, then clustered. It helps defining the amount fo clusters by checking:

% 1. single subject level, which subjects did not have a CAP in a certain
% partition k. It puts a higher boundary of k-selection.
% 2. a notion of how similar areeach subject's CAPs to the concatenated
% average. Can be used as a warning for a subjects.
% 3. each dataset's CAP occurrence rate. It is useful when doing
% test-retest or various sessions with the same subjects.

% The script works for any dimension data, but I suggest you base your %
% analyses using the voxelwise data and clustering results (from either
% percellated, full vw, or downsampled vw).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 19-03-2020
% _________________________________________________________________________

clc,clear
%% 1. EDIT HERE
% path to dFC_CAPs_scripts and external toolboxes
addpath(genpath('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127'))

% 1.a load data, inputs, mask, and parcellation info.
main = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/';
addpath(genpath(main))
main_analysis = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/analysis_vw_plus';

ndatasets = 2; % define the amount of datasets to compare
ncaps = [4:6];    % the clusterings you want to compare... make sure all datasets and runs have them.
analysis_id = 'caps_vw_awk_ane_plus';   % name of the analysis ID you gave in dFC_caps_cluster_mapping_vw_parcel.m

cd(main_analysis)
load('caps_quality')

%Load the results from each datasets for the best run
lead_run = caps_quality.lead_run;
for k = 1:length(ncaps)
    cd(main)
    load([analysis_id '_r' num2str(lead_run(k)) '/' analysis_id '_r' num2str(lead_run(k))])
    results{k} = caps_results.(['caps_' num2str(ncaps(k))]);
    clear caps_results;
end

load([analysis_id '_r' num2str(lead_run(k)) '/inputs'])
ds_name = inputs.ds_name;

% check for removed subjects, and make subject lists
rm_subjects = inputs.rm_subjects;
ss=0;
for nd = 1:ndatasets
    nsubs(nd) = inputs.Nsubs{nd};
    if isempty(rm_subjects{nd}) == 0
        nsubs(nd)=nsubs(nd)-length(rm_subjects{nd});
    end
    for sub = 1:nsubs(nd)
        ss=ss+1;
        sub_list{nd}(sub) = ss;
    end
end

nsubs_concat=ss;

%% 2. Compute within each dataset,k, and cap, the single-subject CAP similarity with the concatenated CAP.
%also compute the mean CAP for each dtaset
cd(main_analysis)
no_cap = cell(length(ncaps),ndatasets);
for k = 1:length(ncaps)
    cmap_concat{k} = results{k}.cap_mean_map;
    cap_ord{k} = caps_quality.c_ind{k}(1,:);
    for nd = 1:ndatasets
        %compute similarity
        no_cap{k,nd}=zeros(nsubs(nd),ncaps(k));
        for c = 1:ncaps(k)
            for sub = 1:nsubs(nd)
                cmap{k,nd}(sub,c,:) = results{k}.cap_mean_map_ss{sub_list{nd}(sub)}(cap_ord{k}(c),:);
                cap_ss_to_concat_sim{k,nd}(sub,c) = corr(spm_vec(cmap_concat{k}(cap_ord{k}(c),:)),spm_vec(cmap{k,nd}(sub,c,:)));
                % ensure all subjects have the cap
                if isnan(sum(cmap{k,nd}(sub,c,:)))
                     no_cap{k,nd}(sub,c) = 1;
                     cmap{k,nd}(sub,c,:)=0;
                     cap_ss_to_concat_sim{k,nd}(sub,c)=0;
                end
                    
            end 
            cap_ss_to_concat_sim_mean{k}(nd,c) = mean(cap_ss_to_concat_sim{k,nd}(:,c));
            cap_ss_to_concat_sim_std{k}(nd,c) = std(cap_ss_to_concat_sim{k,nd}(:,c));
            cap_ds_mean{k}(nd,c,:) = mean(squeeze(cmap{k,nd}(:,c,:)),1);
            for sub = 1:nsubs(nd)
                cap_ss_to_ds_sim{k,nd}(sub,c) = corr(spm_vec(cap_ds_mean{k}(nd,c,:)),spm_vec(cmap{k,nd}(sub,c,:)));
                if isnan(sum(cmap{k,nd}(sub,c,:)))
                     cap_ss_to_ds_sim{k,nd}(sub,c)=0;
                end
            end 
            cap_ss_to_ds_mean{k}(nd,c) = mean(cap_ss_to_ds_sim{k,nd}(:,c));
            cap_ss_to_ds_std{k}(nd,c) = std(cap_ss_to_ds_sim{k,nd}(:,c));            
        end
        no_cap_sum(k,nd) = sum(sum(no_cap{k,nd}));
    end
end

% Plot the within dataset average single-subject to dataset CAP map similarity.
for k = 1:length(ncaps)
    for c = 1:ncaps(k)
        for nd = 1:ndatasets
            y(c,nd) = cap_ss_to_ds_mean{k}(nd,c);
            e(c,nd) = cap_ss_to_ds_std{k}(nd,c)/sqrt(nsubs(nd));
            labelcap{k}{c} = ['C' num2str(c)];
        end
    end
    figure
    fig = gcf;
    fig.Units = 'centimeters';
    fig.Position(3) = 10+k;
    fig.Position(4) = 6;
    fig.Position(1) = 20;
    fig.Position(2) = 20;
    barwitherr(e,y)
    xticks(1:ncaps(k))
    xticklabels(labelcap{k})
    h_legend = legend(ds_name);
    set(h_legend,'Location','bestoutside')

    ylabel({'Average within Dataset';'CAP stability'})
    box off
    fig.PaperPositionMode = 'auto';
    print(['cap_ss_to_ds_sim_k' num2str(ncaps(k))], '-dpng','-r600')
    clear x y e labelcap
    close all
end

%now plot the line with error bar of ALL caps' similarity to the mean
%dataset CAP.
for k = 1:length(ncaps)
        for nd = 1:ndatasets
            y2(nd,k) = mean(cap_ss_to_ds_mean{k}(nd,:));
            e2(nd,k) = std(cap_ss_to_ds_mean{k}(nd,:));
            if isnan(y2(nd,k))
                y2(nd,k)=0;
                e2(nd,k)=0;
            end
            labelk{k} = ['k' num2str(ncaps(k))];
        end
end    
figure
fig = gcf;
for nd = 1:ndatasets
    errorbar(ncaps,spm_vec(y2(nd,:)),spm_vec(e2(nd,:)))
    hold on
end
xlabel('K','fontsize',12)
xlim([0 ncaps(end)+1])
h_legend = legend(ds_name);
set(h_legend,'Location','bestoutside')

ylabel({'Average within Dataset';'CAP stability'})
box off
fig.PaperPositionMode = 'auto';
print(['mean_cap_ss_to_ds_sim_' ds_name{nd} '_k' num2str(ncaps(k))], '-dpng','-r600')
close all

%% 3. compare the CAP occurrence rates
for k = 1:length(ncaps)
    for nd = 1:ndatasets
        for c = 1:ncaps(k)
            for sub = 1:nsubs(nd)
                cap_occ{k}(nd,sub,c) = results{k}.occ_prob_sub(sub,c);
            end
            cap_occ_mean{k}(nd,c)=mean(cap_occ{k}(nd,:,c));
            cap_occ_std{k}(nd,c)=mean(cap_occ{k}(nd,:,c));
            cap_occ_group{k,nd}(:,c) = results{k}.occ_prob_sub(sub_list{nd},c);
            cap_occ_group_mean{k}(nd,c) = mean(cap_occ_group{k,nd}(:,c));
            cap_occ_group_std{k}(nd,c) = std(cap_occ_group{k,nd}(:,c));
        end
    end
    % compare between groups
    for nd1=1:ndatasets
        for nd2=1:ndatasets
            for c = 1:ncaps(k)
                [h_occ{nd1,nd2}{k}(c),p_occ{nd1,nd2}{k}(c),~,~] = ttest2(spm_vec(cap_occ_group{k,nd1}(:,c)),spm_vec(cap_occ_group{k,nd2}(:,c)),'tail','both');
            end
        end
    end
end

%plot CAP-occurrence rates for each k
for k = 1:length(ncaps)
    for c = 1:ncaps(k)
        for nd = 1:ndatasets
            y(c,nd) = cap_occ_group_mean{k}(nd,c);
            e(c,nd) = cap_occ_group_std{k}(nd,c)/sqrt(nsubs(nd));
        end
        labelcap{k}{c} = ['C' num2str(c)];

    end
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 20+k;
    barwitherr(e,y)    
    xticks(1:ncaps(k))
    xticklabels(labelcap{k})
    h_legend = legend(ds_name);
    set(h_legend,'Location','bestoutside')
    ylabel({'Average proportion';'of frames in a CAP'})
    box off
    print(['fprop_cap_k' num2str(ncaps(k))], '-dpng','-r600')
    clear x y e 
    close all
end


%% 4. save data
caps_props.cap_occ_group = cap_occ_group;
caps_props.cap_occ_group_mean = cap_occ_group_mean;
caps_props.cap_occ_group_std = cap_occ_group_std;
caps_props.cap_ord = cap_ord;
caps_props.cap_ss_to_concat_sim = cap_ss_to_concat_sim;
caps_props.cap_ss_to_concat_sim_mean = cap_ss_to_concat_sim_mean;
caps_props.cap_ss_to_concat_sim_std = cap_ss_to_concat_sim_std;
caps_props.cap_ss_to_ds_sim = cap_ss_to_ds_sim;
caps_props.cap_ss_to_ds_sim_mean = cap_ss_to_ds_mean;
caps_props.cap_ss_to_ds_sim_std = cap_ss_to_ds_std;
caps_props.no_cap = no_cap;
caps_props.no_cap_sum = no_cap_sum;
caps_props.cap_occ_group_comparison_pval = p_occ;

save(['caps_props'] ,'caps_props','-v7.3')






end % function;