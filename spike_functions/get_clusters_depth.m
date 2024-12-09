function ids_depth = get_clusters_depth(cluster_info_gen_file)
%% read the file cluster_info.tsv and extract cluster id and depth info
    f=fopen(cluster_info_gen_file);
    %first line is the header
    t=fgetl(f);
    ids_depth=[];
    while ~feof(f)
        t=fgetl(f);
        strline=strsplit(t);
            id=str2num(strline{1});
            depth=str2num(strline{7});
            ids_depth=[ids_depth; [id,depth]];        
    end
    fclose(f);
end