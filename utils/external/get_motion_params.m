cd('D:\Gutierrez\mice_dFC\Xan_on\veh\subjects')
list = dir;
list(1:2) = [];
Nsubs = length(list);

for sub = 1:Nsubs
    sub_list{sub,1} = list(sub).name;
    sub_list1{sub,1} = list(sub).name(1:10);
    sub_list2{sub,1} = list(sub).name(13:25);
end
for sub = 1:Nsubs
    cd(sub_list{sub})
    mkdir('motion_correction')
    cd ..
end

for sub = 1:Nsubs
    cd('D:\Gutierrez\mice_dFC\Xan_on\motion_daniel_scan_10')
    a=dir([sub_list1{sub,1} '*']);
    motion_raw1{sub,1} = load(a.name);
      

    cd('D:\Gutierrez\mice_dFC\Xan_on\motion_daniel_scan11')
    a=dir([sub_list1{sub,1} '*']);
    motion_raw2{sub,1} = load(a.name);

    
    for d = 1:6
        motion_new{sub,1}(1:500,d) = motion_raw1{sub,1}(1001:1500,d);
        motion_new{sub,1}(501:1000,d) = motion_raw2{sub,1}(1:500,d);
    end
    
    cd('D:\Gutierrez\mice_dFC\Xan_on\veh\subjects')
    cd(sub_list{sub})
    cd('motion_correction')
    mot_temp = motion_new{sub};
    name = [sub_list{sub}(1:end-4) '_mcf_02.txt'];
    dlmwrite(name,mot_temp)
end

