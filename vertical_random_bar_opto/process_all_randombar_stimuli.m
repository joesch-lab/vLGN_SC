%script collects all recordings of opto-random vertical bars and starts
%computaion of RFs; one stimulus per recording


%Olga Symonova
%22.09.2020
%
folder_all_bars = 'D:\data\VerticalBars_wOpto';
resroot = 'D:\data\resub\data\vertical_bars_opto';

folder_list=dir(folder_all_bars);

for di=1:length(folder_list)
    if ~isfolder(fullfile(folder_list(di).folder, folder_list(di).name)) || length(folder_list(di).name) < 6 ||...
            isempty(str2num(folder_list(di).name(1:6))), continue; end 
    disp(folder_list(di).name); 
    clearvars -except folder_all_bars resroot folder_list di;
    folder_to_analyze=fullfile(folder_all_bars, folder_list(di).name);

    %collect all needed files automatically, they should be inside folder_to_analyze
    spike_cluster_file= fullfile(folder_to_analyze,'spike_clusters.npy');
    spike_time_file= fullfile(folder_to_analyze,'spike_times.npy');
    cluster_info_file= fullfile(folder_to_analyze,'cluster_group.tsv');
    cluster_info_gen_file= fullfile(folder_to_analyze,'cluster_info.tsv');

    %check if all files are there and return dat and log files
    datfile= check_files_for_analysis(folder_to_analyze, spike_cluster_file, spike_time_file, cluster_info_file);
    if isempty(datfile)
        return;
    end

    %for this session add folder with the functions to read npy files
    addpath(fullfile(pwd,'npy-matlab'));
    
    logdirninfo=dir(fullfile(folder_to_analyze,'stimuli*log.mat'));
    fullfilename_log=fullfile(folder_to_analyze, logdirninfo(1).name);     

    %get the syncing and opto signal       
    [camtrig, red_signal, ~, ~, ~, ~, optotrig, sample_rate]= read_all_from_Intan_RHD2000_file(datfile);
    %create an array with 1 when there was an opto interval    
    [optoHz, opto_pulse_on_s, opto_pulse_cycle_s, opto_cont] = get_continuous_opto(optotrig, sample_rate);
%     figure; plot(optotrig); hold on; plot(opto_cont);
         
    analyze_random_vertical_bar_opto_recording_subsmaple(red_signal, opto_cont, sample_rate, fullfilename_log, resroot);              
end

%% check that all needed files exist in the given folder, issue a warning msg othewise
function datfile = check_files_for_analysis(folder_to_analyze, spike_cluster_file, spike_time_file, cluster_info_file)
%check if the files exist in the folder
datfile='';
logfile='';

if ~exist(spike_cluster_file,'file')
    disp('File spike_clusters.npy is not in the folder. Exitting...');
    return;
end

if ~exist(spike_time_file,'file')
    disp('File spike_times.npy is not in the folder. Exitting...');
    return;
end

if ~exist(cluster_info_file,'file')
    disp('File cluster_info.tsv is not in the folder. Exitting...');
    return;
end

%find dat and log files
d = dir(folder_to_analyze);
nfiles = length(d);

for fi=3:nfiles
    filename=d(fi).name;
    if endsWith(filename,'.dat')
        datfile=fullfile(folder_to_analyze,filename);
    end
    %          if endsWith(filename,'.log')
    %              logfile=fullfile(folder_to_analyze,filename);
    %          end
end
if ~exist(datfile,'file')
    disp('Could not find *.dat file in the folder. Exitting...');
    return;
end
%      if ~exist(logfile,'var')
%         disp('Could not find *.log file in the folder. Exitting...');
%         return;
%      end
end
%%

