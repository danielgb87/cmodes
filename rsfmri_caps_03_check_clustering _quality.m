function caps_03_dFC_caps_clustering_quality

% this function performs various analyses to assess the quality of
% clustering of kmeans given concatenated dataset(s), and various runs of the algorithm
% for each clustering size (k), assessing the algorithm's Variance
% explained, relative distortion. It also computes the between run
% consistency (reproducibility), and how CAPs are resilient to increasing
% cluster size.

% THIS IS A STEP THAT IS PERFORMED ON THE CONCATENATED DATASET, WITHOUT
% TAKING INTO ACCOUNT ANY GROUP-DEPENDENT EFFECTS. IT SERVES AS A PROXY OF
% HOW THE ALGORITHM IS BEING CONSISTENT, AND IS ONLY WORTH IT IF YOU
% PERFORMED VARIOUS RUNS OF CLUSTERING WITH INCREASING K.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Daniel Gutierrez-Barragan 2020, V1 19-03-2020
% _________________________________________________________________________

clc, clear
%% 1. EDIT HERE
% path to dFC_CAPs_scripts and external toolboxes
addpath(genpath('/home/dgutierrez/scripts_toolboxes/Analysis_scripts/caps_scripts_210127'))

% 1.a load results, inputs, other info from folder containing the directories of clustering results.
main = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221';
addpath(genpath(main))
main_analysis = '/media/DATA/dgutierrez/CAPS_tutorial/caps_scripts_concat_230221/analysis_vw_plus';

cd(main)
mkdir('analysis_vw_plus')

ndatasets = 2; % define the amount of datasets to compare
ncaps = [4:6];    % the clusterings you want to compare... make sure all datasets and runs have them.
is_lin_inc = 1; % state if the amount of clusters increases in k=1 units (ex. ncaps=2:1:12)

%1.b define the following for each run (in general, all runs would have the same ID, with a termination _r1, run number):
analysis_id = 'caps_vw_awk_ane_plus';   % name of the analysis ID you gave in dFC_caps_cluster_mapping_vw_parcel.m
analysis_name = 'awk_ane_scrub';
nruns = 2;    % number of runs 

%Load the results from each dataset and run
for r = 1:nruns
        cd(main)
        load([analysis_id '_r' num2str(r) '/' analysis_id '_r' num2str(r)])
        results{r} = caps_results;
        clear caps_results;
end



%% 2. Compute within a dataset: Variance explained and Relative distortion curves
cd(main_analysis)
for k = 1:length(ncaps)
    for r = 1:nruns
        [Iw(r,k),Ib(r,k),VE(r,k),RD(r,k)] = dFC_utils_caps_variance_explained(results{r}.(['caps_' num2str(ncaps(k))]));
    end
    
    VE_mean(k) = mean(VE(:,k));
    RD_mean(k) = mean(RD(:,k));
    VE_std(k) = std(VE(:,k));
    RD_std(k) = std(RD(:,k));
    [~, lead_run(k)] = max(VE(:,k));

    display(['done with k=' num2str(ncaps(k))])
end

% plot results
% variance explained and relative distortion 
figure
fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 15; fig.Position(4) = 6;
fig.Position(1) = 5; fig.Position(2) = 5;
shadedErrorBar(ncaps,VE_mean(:),VE_std(:),'b')
ylabel('Average variance explained')
xticks(1:ncaps(end))
xlabel('Number of clusters')
xlim([0 ncaps(end)+1])
fig.PaperPositionMode = 'auto';
print(['VE_' analysis_name], '-dpng','-r600')

figure
fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 15; fig.Position(4) = 6;
fig.Position(1) = 5; fig.Position(2) = 5;
shadedErrorBar(ncaps,RD_mean(:),RD_std(:))
ylabel('Average relative distortion')
xticks(1:ncaps(end))
xlabel('Number of clusters')
xlim([0 ncaps(end)+1])
fig.PaperPositionMode = 'auto';
print(['RD_' analysis_name], '-dpng','-r600')
close all

% 2.2 Between run stability
cap_sim = cell(length(ncaps),1);
c_ind = cell(length(ncaps),1);

for k = 1:length(ncaps)
    % match each CAP from each run to the one in the Leading run.
    tmp = cell(nruns,1);
    [~, run_ind(k,:)] = sort(squeeze(VE(:,k)),'descend');
    for r = 1:nruns
        tmp{r} = results{run_ind(k,r)}.(['caps_' num2str(ncaps(k))]);
    end
    Dmat = cell(nruns);
    [~, c_ind{k}(1,:)] = sort(tmp{1}.occ_prob_mean,'descend');
    
    for r1 = 2:nruns
        for i = 1:ncaps(k)
            for j = 1:ncaps(k)
                Dmat{r1}(i,j) = pdist([(spm_vec(tmp{1}.cap_mean_map(c_ind{k}(1,i),:))'); (spm_vec(tmp{r1}.cap_mean_map(j,:))')]);
            end
        end
        [c_ind{k}(r1,:), ~] = munkres_HA(Dmat{r1});
    end
    % compare caps from runs
    for c = 1:ncaps(k)
        count = 1;
        for r1 = 1:nruns-1
            for r2 = r1+1:nruns              
                cap_sim{k}(count,c) = corr(tmp{r1}.cap_mean_map(c_ind{k}(r1,c),:)', tmp{r2}.cap_mean_map(c_ind{k}(r2,c),:)');
                count=count+1;
            end
        end
        cap_sim_mean{k}(c) = mean(cap_sim{k}(:,c));
        cap_sim_std{k}(c) = std(cap_sim{k}(:,c));
    end   
    clear Dmat tmp    
end

% plot between run stability
for k = 1:length(ncaps)
    % plot and save cap similarity in the mixed data analysis.
    for c = 1:ncaps(k)
        y(c,1) = cap_sim_mean{k}(c);
        e(c,1) = cap_sim_std{k}(c);
        labelcap{c} = ['C' num2str(c)];
    end
    figure
    fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 10+ncaps(k);
    fig.Position(4) = 6; fig.Position(1) = 5; fig.Position(2) = 5;
    barwitherr(e,y)
    hold on
    hline = refline([0 0.8]);
    hline.Color = 'r';
    set(gca,'XTickLabel',labelcap)
    xticks(1:ncaps(k))
    ylabel('mean between run stability')
    ylim([-0.2 1.2])
    box off
    fig.PaperPositionMode = 'auto';
    print(['btw_run_stability_' analysis_name '_k' num2str(ncaps(k))], '-dpng','-r600')
    clear x y e labelcap
    close all
end

% 2.3 Compute within each run, the evolution of CAPs as k grows. this is
% only posible if there are more than 2 k's used, and if they are
% increasing by 1 cluster

if is_lin_inc==1
    cap_corr = cell(nruns,1);
    cap_evo = cell(nruns,1);
    cap_ord = cell(nruns,1);
    for r = 1:nruns
        cap_corr{r} = cell(length(ncaps),1);
        cap_evo{r} = cell(length(ncaps),1);
        cap_ord{r} = cell(length(ncaps),1);
        cap_corr{r}{1} =zeros(1,min(ncaps));
        cap_ord{r}{1} = c_ind{1}(r,:);
        for k = 2:length(ncaps)
            tmp1 = results{r}.(['caps_' num2str(ncaps(k-1))]).cap_mean_map;
            tmp2 = results{r}.(['caps_' num2str(ncaps(k))]).cap_mean_map;
            % assign the caps from k to the k-1 clustering
            for i = 1:ncaps(k)-1
                for j = 1:ncaps(k)
                    Dmate{r}{k}(i,j) = pdist([spm_vec(tmp1(cap_ord{r}{k-1}(i),:))'; spm_vec(tmp2(j,:))']);
                end
            end
            [cap_ord{r}{k}, ~] = munkres_HA(Dmate{r}{k});
            % identify the "new CAP" in k
            v = 1:ncaps(k);
            for c = 1:ncaps(k-1)
                if ismember(cap_ord{r}{k},v(c)) ==0
                    cap_ord{r}{k}(ncaps(k-1)+1) = c;
                end
            end

            %compute the similarities between matched caps
            for c = 1:ncaps(k-1)
                cap_corr{r}{k}(c) = corr(spm_vec(tmp1(cap_ord{r}{k-1}(c),:)), spm_vec(tmp2(cap_ord{r}{k}(c),:)));
            end
            cap_corr{r}{k}(end+1) = 0;

        end
        % organize evolution vectors
        for k = 1:length(ncaps)
            for kk = k:length(ncaps)
                cap_evo{r}{k}(kk) = cap_corr{r}{kk}(k);
            end
        end
        figure
        fig = gcf; fig.Units = 'centimeters'; fig.Position(3) = 35;
        fig.Position(4) = 10; fig.Position(1) = 0; fig.Position(2) = 0;
        plot(2:length(ncaps)+1,cap_evo{r}{1}(1:length(ncaps)));
        hold on
        for k = 2:length(cap_evo{r})
            plot(k:length(ncaps)+1,cap_evo{r}{k}(k-1:length(ncaps)))
            hold on
        end
        hold on
        xlabel('K','fontsize',12)
        ylabel('CAP stability','fontsize',12)
        xlim([0 ncaps(end)+1])
        ylim([0 1])
        fig.PaperPositionMode = 'auto';
        print(['wtn_cap_stability_evolution_' analysis_name '_run' num2str(r)], '-dpng','-r600');
        close all
    end
end

%% 3. organize results
caps_quality.lead_run = lead_run;
caps_quality.c_ind = c_ind;
if is_lin_inc==1
    caps_quality.cap_corr = cap_corr;
    caps_quality.cap_evo = cap_evo;
end
caps_quality.cap_sim = cap_sim;
caps_quality.cap_ord = cap_ord;
caps_quality.Ib=Ib;
caps_quality.Iw=Iw;
caps_quality.VE=VE;
caps_quality.RD=RD;
caps_quality.run_ind=run_ind;

save(['caps_quality'] ,'caps_quality','-v7.3')
    
       




end % function