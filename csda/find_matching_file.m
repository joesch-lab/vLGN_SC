function csdfile = find_matching_file(folder2searh, filename)
%find csd-file
[~,recname,~]=fileparts(filename);
recwords=strsplit(recname,{'_','-'});
numr=cellfun(@numel, recwords);
recwords=recwords(numr>1);

%find best matching csd file
maxmatch=0;
matid=-1;
csda_list=dir(fullfile(folder2searh,'*.mat'));
for ci=1:length(csda_list)
    maxmatchi=0;
    for wi=1:length(recwords)
        if contains(csda_list(ci).name,recwords(wi))
            maxmatchi=maxmatchi+1;
        end
    end
    if maxmatchi>=maxmatch
        maxmatch=maxmatchi;
        matid=ci;
    end
end

csdfile=fullfile(csda_list(matid).folder, csda_list(matid).name);
disp(['Searching csd info for ', filename]);
disp(['matching file found: ', csdfile]);
end