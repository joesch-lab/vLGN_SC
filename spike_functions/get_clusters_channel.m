function ids_cha = get_clusters_channel(cluster_info_gen_file)
%% read the file cluster_info.tsv and extract cluster id and depth info
    f=fopen(cluster_info_gen_file);
    %first line is the header
    t=fgetl(f);
    ids_cha=[];
    while ~feof(f)
        t=fgetl(f);
        strline=strsplit(t);
            id=str2num(strline{1});
            cha=str2num(strline{6});
            ids_cha=[ids_cha; [id,cha]];        
    end
    fclose(f);
end