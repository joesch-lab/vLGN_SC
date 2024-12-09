function log_file_info = parse_log_files(folder_to_analyze, red_signal, sr,start_experiment_s)
% uses the red syncing signal and info stored in stimuli files to parse one
% recording into stimuli.
% params:
%folder_to_analyze: where all the recoded data and stim files are
% red_signal: syncing array of 1 and 0, loaded previously
% sr: sampling rate, loaded together with the triggers
% start_experiment_s: start of the experiment in secs if different from 0

if ~exist('start_experiment_s','var') || isempty(start_experiment_s)
    start_experiment_s=0;
end

log_file_info={};
log_list=dir(fullfile(folder_to_analyze, 'stimuli_*'));

%sort log files by date
[~,idx] = sort([log_list.datenum]);
log_list=log_list(idx);

nfiles=length(log_list);
start_experiment_sample=max(1,start_experiment_s*sr);
for fi=1:nfiles
    stimlogfile=fullfile(log_list(fi).folder, log_list(fi).name);
    log_file_info(fi).name=log_list(fi).name;
    log_file_info(fi).folder=log_list(fi).folder;
    log_file_info(fi).file_tstamp=log_list(fi).date;    
    log_file_info(fi).params = read_stimulus_file(stimlogfile);    
    [start_i, end_i] = start_end_stimulus(red_signal(start_experiment_sample:end),fi,log_list, 0, sr);    
    log_file_info(fi).start_s=double(start_i+start_experiment_sample-1)/sr;
    log_file_info(fi).end_s=double(end_i+start_experiment_sample-1)/sr;
end
%pause(10);
close all;