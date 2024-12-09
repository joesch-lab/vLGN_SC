% this script will read all basic info about experiment and save it into a
% structure to speed up further proccessing

% for a folder with the experiment,
%  collect opto, visual, loco, camera and spiking info

%the txt file list_process lists the names of the folders in the repo to
%process. If such file is not provided all folders in the the repo will be
%processed.

%laser_yn=1 idicates whether optogenetics was used in the experiment

%analyze_ephys=1 if spiking info is available; analyze_ephys=0 if only behav

%Olga Symonova
%17.01.2023
%
% repo_folder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\BackupRecordings\Quilo\KilosortingVZ\Sinusoidal_OKR';
% repo_folder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\BackupRecordings\Quilo\KilosortingVZ\Sinusoidal_OKR';
% folder_proc =  'D:\data\tnt_data';
% csda_resfolder = 'D:\data\tnt_data\csda';

% repo_folder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\Opto_Behavior\IGL';
% folder_proc =  '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\Opto_Behavior\IGL\proc';
% list_process = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\Opto_Behavior\IGL\all_rec_list.txt';
% laser_yn=1;

% list_process = 'D:\data\NewMap\all_rec_list.txt';
% list_process = 'D:\NewMap\with_laser_list.txt'; %laser_yn=1;
% list_process='D:\NewMap\control_list.txt'; % laser_yn=0;
% list_process=[];
% laser_yn=0;

% repo_folder = 'D:\VerticalBars_wOpto';
% folder_proc =  'D:\VerticalBars_wOpto\res';
% csda_resfolder = 'D:\VerticalBars_wOpto\csda_res';
% % list_process='D:\VerticalBars_wOpto\no_laser_list.txt';
% list_process='D:\VerticalBars_wOpto\with_laser_list.txt';
% laser_yn=1;

repo_folder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\Opto_Behavior\in_SC\Inhibitory_Opto';
folder_proc =  '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\Opto_Behavior\in_SC\Inhibitory_Opto\proc_1sec';
list_process = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\Opto_Behavior\in_SC\Inhibitory_Opto\all_rec_list_1sec.txt';
laser_yn=1;

analyze_ephys=1;
 
folderlist = get_list_from_file(list_process);
% 
% if isempty(list_process)
%     folderlist={};
%     list_process=dir(repo_folder);
%     laser_yn=0;
%     for fi=1:length(list_process)
%         if startsWith(list_process(fi).name,'.')
%             continue;
%         else
%             folderlist=[folderlist;list_process(fi).name];
%         end
%     end
% end
% 
% 
% folderlist={'Sham4_MJ1F23-720_1400um_OlgaList'};
% folderlist={'Sham4_MJ1F23-720_1400um_SinCheckers1'; 'Sham4_MJ1F23-720_1400um_SinCheckers2'};

nf=length(folderlist);


%% loop to process new files
for fi=1:nf
    expfolder=folderlist{fi}
    folder_to_analyze =fullfile(repo_folder,expfolder);
    folder_to_save =  fullfile(folder_proc,expfolder);
    mkdir(folder_to_save);
    save_file=fullfile(folder_to_save ,'data.mat');

    data={};

    %% get data from the folder name
    [repopath,fname,~]=fileparts(folder_to_analyze);
    data.exp_info = extract_data_from_foldername(folder_to_analyze);
    data.exp_info.fullpath =folder_to_analyze;%repopath;
    data.exp_info.laser_yn=laser_yn;
    
    %% get data from the rhd file - assume only 1 is in the folder
    datpatt=fullfile(folder_to_analyze,'*.dat');
    datinfo=dir(datpatt);
    if isempty(datinfo)
        datpatt=fullfile(folder_to_analyze,'*.rhd');
        datinfo=dir(datpatt);
    end
    if length(datinfo)~=1
        warning('More than 1 dat file in the folder. Stop and think.');
        continue;
    end
    data.exp_info.dat_file = datinfo(1).name;
    datfile=fullfile(folder_to_analyze,datinfo(1).name);
    data.trigger_data = extract_dat_file_info(datfile);

    %% parse recording among the stimuli files, identifying the start and end of each stimulus by hand,
    % load the params of each stimulus
    data.exp_info.logfiles  = parse_log_files(folder_to_analyze, data.trigger_data.red_signal, data.trigger_data.sampling_rate, 1);
    
    %% get spike sorting data - one set per folder
    if analyze_ephys
        spike_info = get_spike_info(folder_to_analyze);
        data.spike_data=spike_info;
    end
    
    %% if bin file is available, do current density analysis
    bin_file = find_bin_file(folder_to_analyze);
    data.trigger_data.binfile=bin_file;
    if ~isempty(bin_file) %% and opto flash stimulus
        %find id of the local_flash stimuli and use it for csda
        logid = find(contains({data.exp_info.logfiles(:).name},'local_flash'));        
        if ~isempty(logid)
            ststim=uint64(data.exp_info.logfiles(logid).start_s*data.trigger_data.sampling_rate);
            enstim=uint64(data.exp_info.logfiles(logid).end_s*data.trigger_data.sampling_rate);        
            data.csda = CSDanalysis(data, csda_resfolder,1,1,ststim,enstim);
            data = set_csda_depth_per_cluster(data);
        else %find the csda file
            csda_file = find_matching_file(csda_resfolder, data.exp_info.fullpath);
            load(csda_file);
            data.csda=csda_data;
            data = set_csda_depth_per_cluster(data);
        end
    end

    %% DeepLabCut data on the eye
    eye_file = find_dlc_eye_file(folder_to_analyze);
    if ~isempty(eye_file)
        pupil_el_az_area = get_eye_values_from_dlc_file(eye_file);
        data.dlc.pupil_el = pupil_el_az_area(:,1);
        data.dlc.pupil_az = pupil_el_az_area(:,2);
        data.dlc.pupil_area = pupil_el_az_area(:,3);
    end    
    
    
    %% get motion energy using SVD decomposition of difference of frames
    motion_data = get_motion_energy(data, folder_to_save);
    data.motion_energy=motion_data;

    %% (re)save data
    save(save_file,"data","-v7.3");
end




















