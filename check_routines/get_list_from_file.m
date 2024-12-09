function folderlist = get_list_from_file(folderlist_file)
% get the list of files to process saved as rows in the txt file

folderlist ={};
f=fopen(folderlist_file);
i=1;
tline=fgetl(f);
while  ~feof(f)
    if ~isempty(tline)
        folderlist{i} =tline;
        tline=fgetl(f);
        i=i+1;
    end
end
if ~isempty(tline)
    folderlist{i} =tline;
end
end