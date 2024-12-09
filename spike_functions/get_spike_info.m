function spike_info = get_spike_info(exp_folder)
    %load the data after running kilosort
    spike_info={};

    %collect all spike info files automatically
    filelist = dir(fullfile(exp_folder, '*.*'));  %get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);     
    
    id = find(contains({filelist(:).name},'spike_times.npy'));
    if isempty(id)
        warning([exp_folder, ' has no spike data.']);
        return;
    end

    spike_time_file=filelist(id);
    spike_time_file=fullfile(spike_time_file.folder,spike_time_file.name);

    id = find(contains({filelist(:).name},'spike_clusters.npy'));    
    spike_cluster_file=filelist(id);
    spike_cluster_file=fullfile(spike_cluster_file.folder,spike_cluster_file.name);


    id = find(contains({filelist(:).name},'cluster_group.tsv'));
    cluster_info_file=filelist(id);
    cluster_info_file=fullfile(cluster_info_file.folder,cluster_info_file.name);

    id = find(contains({filelist(:).name},'cluster_info.tsv'));
    cluster_info_gen_file=filelist(id);
    cluster_info_gen_file=fullfile(cluster_info_gen_file.folder,cluster_info_gen_file.name);
        
    %get depth info
    ids_depth = get_clusters_depth(cluster_info_gen_file);

    %get channel info
    ids_cha = get_clusters_channel(cluster_info_gen_file);
    
    %cluster label and time of each spike
    cluster_id = readNPY(spike_cluster_file);
    spike_t = readNPY(spike_time_file);

    %find ids of good clusters
    goodclusters = get_clusters_by_property(cluster_info_file,'good');
    muaclusters = get_clusters_by_property(cluster_info_file,'mua');    
    
    n_gcl=length(ids_depth);
    i=1;
    spike_info={};
    for ci=1:n_gcl
        %cluster id as after kilosort
        idc = ids_depth(ci,1);
        if ismember(idc,goodclusters)
            classification = 'good';
        elseif ismember(idc,muaclusters)
            classification = 'mua';
        else
            continue;
        end        
        %select spikes per cluster
        spikes=find(cluster_id==idc);
        %get the time of spikes
        ts=spike_t(spikes); 
        %store results in the structure
        spike_info(i).id=idc;
        spike_info(i).depth=ids_depth(ci,2);         
        spike_info(i).channel = ids_cha(ci,2);
        spike_info(i).class=classification;
        spike_info(i).nspikes=length(spikes);
        spike_info(i).spike_timing=ts;        
        i=i+1;
    end
end