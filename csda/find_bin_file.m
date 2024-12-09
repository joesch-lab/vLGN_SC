function binfile = find_bin_file(expfolder)
    binfile=[];
    bin_pattern = '*.bin';
    bin_mask = fullfile(expfolder,bin_pattern);
    bin_file_info = dir(bin_mask);
    if ~isempty(bin_file_info)
        binfile = fullfile(bin_file_info(1).folder,bin_file_info(1).name);
    end
end
