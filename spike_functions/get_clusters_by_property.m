function ids = get_clusters_by_property(filename, propstr)
%% read the file that encodes the cluster id and associated property
%  returns the indices of the clusters with the given property 'propstr'
    f=fopen(filename);
    %first line is the header
    t=fgetl(f);
    ids=[];
    while ~feof(f)
        t=fgetl(f);
        strline=strsplit(t);
        if strcmp(strline{2},propstr)
            id=str2num(strline{1});
            ids=[ids, id];
        end
    end
    fclose(f);
end
