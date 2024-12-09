% loads all the files with laser and computes responses to each opto-flash
% trial

%Olga Symonova
%14.02.2023

folder_proc =  'D:\data\resub\data';
% list_to_process='D:\NewMap\with_laser_list.txt';
list_to_process = 'D:\data\resub\data\all_rec_list.txt';
% list_to_process='D:\NewMap\control_list.txt';
folderlist = get_list_from_file(list_to_process);
plot_histogram_yn=0;

nf=length(folderlist);
%% loop for pre-processed files
for fi=1:1%nf    
 expfolder=folderlist{fi} 
 folder_to_save =  fullfile(folder_proc,expfolder);
 save_file=fullfile(folder_to_save ,'data.mat');
 load(save_file); %should have data sturcture with trigger data now

 organize_spikes_by_offsets(data, plot_histogram_yn, folder_to_save);
 ofdata_file = fullfile(folder_to_save,'opto_flash_data.mat');
 compute_response_params(ofdata_file, folder_to_save);
end



















