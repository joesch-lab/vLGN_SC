%this script will collect all the files with CSD analysis in the list given
%in the txt file list_process
% it will center all recondings around the inflection point and compute the
% average envelope indicating the extend of visual/opto stimulation.


repo_folder =  'D:\data\resub\data';
fig_folder = 'D:\data\resub\figs';


list_process='D:\data\resub\data\with_laser_list.txt';
folderlist = get_list_from_file(list_process);

nf=length(folderlist);
sr=20000;
probe_all=[];
rec_list={};
%% loop for pre-processed files
for fi=1:nf    
    expfolder=folderlist{fi}    
    folder_to_save =  fullfile(repo_folder,expfolder);
    save_file=fullfile(folder_to_save ,'data.mat');
    load(save_file);
    rec_list(fi)={expfolder};

    probe_i=NaN(64,2);

    if data.exp_info.laser_yn==0, continue; end
    ncha=size(data.csda.CSD_vis_on,1);
    nt=size(data.csda.CSD_vis_on,2);
    onset=data.csda.flash_onset_samples;
    
    cha_inf_pt = find(data.csda.channel_labels==data.csda.inflection_cha_above)+1;

    var_vis=data.csda.envelope_vis./max(data.csda.envelope_vis);
    var_opto=data.csda.envelope_opto./max(data.csda.envelope_opto);
    
    st_cha= (32-(cha_inf_pt-1))+1;    
    en_cha=32+(32-cha_inf_pt)+1;
    probe_i(st_cha:en_cha,1)=var_vis;
    probe_i(st_cha:en_cha,2)=var_opto;
    probe_all=cat(3,probe_all,probe_i);
end

probe_mean=mean(probe_all,3,"omitnan");
%set to 0 channels with no input for plotting
nancha=isnan(probe_mean(:,1));
probe_mean(nancha,1)=0;
nancha=isnan(probe_mean(:,2));
probe_mean(nancha,2)=0;

fig =figure; patch(-probe_mean(:,1), -32:31,[204/255,1,204/255]);
all_depth_um = -32*25:25:31*25;
hold on; patch(probe_mean(:,2),-32:31,[204/255,204/255,1]);
set(gca,'YDir','reverse');
maxcount=max(probe_mean(:));
xlim([-1;1]);
yticks(-32:8:31);
yticklabels(-32*25:8*25:31*25);
hold on; plot([-maxcount-5;maxcount+5],[0,0],'--','Color',[0.2 0.2 0.2])
title('mean csd variance on vis & opto stimulation');
filename=fullfile(fig_folder,'csd_var_envelopes.svg');
saveas(fig,filename);

data_file=fullfile(fig_folder,'csd_var_envelopes_data.mat');
save(data_file,'probe_all','probe_mean','rec_list');


