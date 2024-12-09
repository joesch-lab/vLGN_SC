function eye_dlc = find_dlc_eye_file(expfolder)
    eye_dlc=[];
    eye_dlc_pattern = '\**\*eye.csv';
    eye_file_mask = fullfile(expfolder,eye_dlc_pattern);
    eye_file_info = dir(eye_file_mask);
    if ~isempty(eye_file_info)
        eye_dlc = fullfile(eye_file_info(1).folder,eye_file_info(1).name);
    end
end

