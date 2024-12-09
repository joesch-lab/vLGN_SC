function datafilelist = collect_processed_stimuli_files(resfolder,stimname, musthave_dlc_info, musthave_spike_info, musthave_loco_info)
% crawls thru resfolder and collects all *data.mat files with given
% stimulus name that have various requested information

datafilelist={};
k=1;
filelist=dir(fullfile(resfolder,['**',filesep,'*data.mat']));
for fi=1:length(filelist)
    clear data;
    datafile=fullfile(filelist(fi).folder, filelist(fi).name);
    load(datafile);

    %skip recordings the stimuli
    if ~exist('data',"var") ...
        || ~(isfield(data, 'exp_info') && isfield(data.exp_info, 'logfiles') &&...
        ~isempty(data.exp_info.logfiles) && contains([data.exp_info.logfiles(:).name],stimname)) 
        continue;
    end

    if musthave_dlc_info && (~isfield(data, 'dlc') || isempty(data.dlc) )
        disp(['File ',datafile, 'has no dlc info']);
        continue;
    end

    if musthave_spike_info && (~isfield(data, 'spike_data') || isempty(data.spike_data))
        disp(['File ',datafile, 'has no spike info']);
        continue;
    end

    if musthave_loco_info && (~isfield(data.trigger_data, 'loco') || isempty(data.trigger_data.loco))
        disp(['File ',datafile, 'has no loco info']);
        continue;
    end
    
    datafilelist{k}=datafile;
    k=k+1;
end
end