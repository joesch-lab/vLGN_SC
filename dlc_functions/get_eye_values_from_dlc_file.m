function pupil_el_az_area = get_eye_values_from_dlc_file(dlc_eye_file)

    %first clean the file to remove comas in the names, otherwise reading
    %the file will produce errors
    dlc_file_new = fix_file_names(dlc_eye_file);    
    opts = detectImportOptions(dlc_file_new);
    opts.SelectedVariableNames = {'pupil_width','pupil_height','pupil_el','pupil_az'};
    pupil_vals = readmatrix(dlc_file_new, opts);
    pupil_area = pi*pupil_vals(:,1).*pupil_vals(:,2)/4;
    pupil_el_az_area = cat(2,pupil_vals(:,3:4),pupil_area); 
end