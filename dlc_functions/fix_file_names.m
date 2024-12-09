function dlc_file_new = fix_file_names(dlc_file)
%some file names contain comas, then the field inside csv file with labels
%will contain comas, then it will be read incorrectly, as it is assumed
%that comas are used to separate values.

%before importing the file, check if the value tiff_name which is contained
%between | | has commas. if so remove them


dlc_file_new=[dlc_file(1:end-4),'_clean.csv'];
if exist(dlc_file_new,"file")~=2
    fid =fopen(dlc_file);
    fid_new =fopen(dlc_file_new,'w');
    %first three lines are headers
    tline = fgetl(fid);
    while ischar(tline)
        vdelim_pos = strfind(tline,'|');
        if ~isempty(vdelim_pos) && length(vdelim_pos>1)
            subline = tline(vdelim_pos(1):vdelim_pos(2));
            subline = strrep(subline, ',','_');
            tline(vdelim_pos(1):vdelim_pos(2))=subline;
        end
        fprintf(fid_new, tline);
        fprintf(fid_new, '\n');
        tline = fgetl(fid);
    end
    fclose(fid);
    fclose(fid_new);
end
end