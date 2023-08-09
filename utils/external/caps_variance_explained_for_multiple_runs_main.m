% in this part the script calculates the caps variance explained for each
% runs

%here you need to insert the path to the folder with the different runs
base_path = '/media/DATA1/Fabio_Gatto_studies/MRC_Human_rsfMRI/CAPs_analysis_pp/CAPs_with_parcelletions/MRC_Craddock_10_runs/';
cd(base_path)
%insert the number of runs
n_runs = 10;
VE_TD_runs=cell(n_runs,1);
VE_AD_runs=cell(n_runs,1);
%the following is the loop for TD subjects
for i=1:n_runs
    results_folder = ['MRC_TD_Craddock_950parcel_r' num2str(i) '_plus_k2-40'];
    cd(results_folder)
    load('caps_results.mat');
    N_caps=length(caps_results);
    Ib_vec=NaN([N_caps,1]);
    Iw_vec=NaN([N_caps,1]);
    VE_vec=NaN([N_caps,1]);
    RD_vec=NaN([N_caps,1]);
    for kk=1:N_caps
        results=caps_results{1,kk};
        [Iw,Ib,VE,RD] = caps_utils_compute_varexp(results);
        Iw_vec(kk,1)=Iw;
        Ib_vec(kk,1)=Ib;
        VE_vec(kk,1)=VE;
        RD_vec(kk,1)=RD;
    end
    filename=sprintf('VE_TD_Craddock_parcel_r%d.mat',i);
    save(filename,'VE_vec','-v7.3');
    VE_TD_runs{i}=VE_vec;
    cd ..
end

%now there is the loop for the AD population
for r=1:n_runs
    results_folder = ['MRC_AD_Craddock_950parcel_r' num2str(r) '_plus_k2-40'];
    cd(results_folder)
    load('caps_results.mat');
    N_caps=length(caps_results);
    Ib_vec=NaN([N_caps,1]);
    Iw_vec=NaN([N_caps,1]);
    VE_vec=NaN([N_caps,1]);
    RD_vec=NaN([N_caps,1]);
    for rr=1:N_caps
        results=caps_results{1,rr};
        [Iw,Ib,VE,RD] = caps_utils_compute_varexp(results);
        Iw_vec(rr,1)=Iw;
        Ib_vec(rr,1)=Ib;
        VE_vec(rr,1)=VE;
        RD_vec(rr,1)=RD;
    end
    filename=sprintf('VE_AD_Craddock_parcel_r%d.mat',r);
    save(filename,'VE_vec','-v7.3');
    VE_AD_runs{r}=VE_vec;
    cd ..
end
%here you save the two cell array with the VE for each run
save('VE_AD_runs.mat','VE_AD_runs','-v7.3');
save('VE_TD_runs.mat','VE_TD_runs','-v7.3');

k=2:1:40;
k=k';
M_AD=zeros(39,1);
for i=1:10
    M_AD (:,1)= M_AD(:,1) + VE_AD_runs{i}(:,1);
end
M_AD = M_AD/10;

AD=NaN(39,10);
for i=1:10
    AD(:,i) = VE_AD_runs{i}(:,1);
end

TD=NaN(39,10);
for i=1:10
    TD(:,i) = VE_TD_runs{i}(:,1);
end



S_AD=std(AD, 0, 2);
M_AD=mean(AD,2);
S_TD=std(TD, 0, 2);
M_TD=mean(TD,2);

f1=figure;
shadedErrorBar(k,M_TD,S_TD,'b');
ylabel('VE_TD');
xlabel('k');
title('Variance Explained of TD subjects for 10 runs');
f2=figure;
shadedErrorBar(k,M_AD,S_AD,'r');
ylabel('VE_AD');
xlabel('k');
title('Variance Explained of AD subjects for 10 runs');
f3=figure;
shadedErrorBar(k,M_TD,S_TD,'b');
ylabel('VE');
xlabel('k');
% legend('TD','AD');
title('Variance Explained of TD and AD subjects for 10 runs');
hold on
shadedErrorBar(k,M_AD,S_AD,'r');
%hold on 
%plot([12 12], [0.58 0.6855],'b','LineWidth',1.5);
print(f3,'VE_10_runs.png', '-dpng', '-r300');
% in this second part you calculate the highest curvature for TD subjects
k_TD_curv=k(2:15);
M_TD_curv=M_TD(2:15)';
Vertices_TD=[k_TD_curv M_TD_curv];
Curv_TD=LineCurvature2D(Vertices_TD);
Curvature_TD=[k_TD_curv Curv_TD];
f4=figure;
plot(Curvature_TD(:,1),Curvature_TD(:,2),'b','LineWidth',1.2);
ylabel('Curvature_TD(k)');
xlabel('k');
title('Maximum Curvature of TD subjects');
%hold on
%plot([12 12], [Curv_TD(1) Curv_TD(10)],'r','LineWidth',1.5);
print(f4,'Maximum_Curvature_TD.png', '-dpng', '-r300');
% in this second part you calculate the highest curvature for AD subjects
k_AD_curv=k(2:15);
M_AD_curv=M_AD(2:15)';
Vertices_AD=[k_AD_curv M_AD_curv];
Curv_AD=LineCurvature2D(Vertices_AD);
Curvature_AD=[k_AD_curv Curv_AD];
f5=figure;
plot(Curvature_AD(:,1),Curvature_AD(:,2),'b','LineWidth',1.2);
ylabel('Curvature_AD(k)');
xlabel('k');
title('Maximum Curvature of AD subjects');
%hold on
%plot([10 10], [-0.004 Curv_AD(8)],'r','LineWidth',1.5);
print(f5,'Maximum_Curvature_AD.png', '-dpng', '-r300');

% alternative method to calculate the elbow of the curve
[k_max_curv_TD, idx_of_result_TD] = knee_pt(M_TD,k);
[k_max_curv_AD, idx_of_result_AD] = knee_pt(M_AD,k);



f6=figure;
plot(k,M_TD,'b');
ylabel('VE');
xlabel('k');
title('Variance Explained of TD and AD subjects for 10 runs');
hold on
plot(k,M_AD,'r');
legend('TD','ASD');
%hold on 
%plot([12 12], [0.58 0.6855],'b','LineWidth',1.5);
print(f6,'VE_10_runs_normal_plot.png', '-dpng', '-r300');