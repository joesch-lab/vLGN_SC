%% bahavioral analysis

% %from all the recordings with opto-flashes select those of the given length, laser yes/no
% collect responses of all animals and reps

% select opto flashes (maximally) separated from other opto+visual flashes
% select the time before opto flash
% get the loco and pupil, face energy for the time before+opto
% separate all reps into those with stationary animal, vs. moving
% plot: behav traces forward, sideways, rotation, for stationary and moving
% plot: scatter pupil R vs az, pupil R vs el, diff colors for before and
% opto, check scatter for before/opto+still, before/opto+moving
% plot: facial components before/opto for still &moving
% 

path_to_figures = 'D:\data\resub\figs';

control=1;
flash_dur_sec=1; %0.2; %0.2 %1;
file_str = 'optoloco_';
if control
    file_str=[file_str,'control'];
end
if flash_dur_sec==1
    file_str=[file_str,'_tflash_1s'];
else
    file_str=[file_str,'_tflash_200ms'];
end

if control
    repo_folder='D:\data\resub\data\control';
    folderlist={'WT_TaperedFiber_MJ1F23-874_LocalFlash_Opto_LPw10_1',...
        'WT_TaperedFiber_MJ1F23-875_LocalFlash_Opto_LPw10_1', 'WT_TapFiber_vLGN_M1F23-874_LocalFlash_Opto_1sec_LPwr5_1'};
else
    repo_folder='D:\data\resub\data\';
    list_process='D:\data\resub\data\with_laser_list.txt';
    folderlist = get_list_from_file(list_process);
end


nf=length(folderlist);
sr = 20000;

%% set params for processing based on the flash duration
if flash_dur_sec==1
    %% values for the analysis of 1sec flash
    flash_dur=flash_dur_sec*sr;
    tbefore=flash_dur;
    tafter=tbefore;

    %azimuth vals analysis
    t_start_before_azi = 0.25*sr; %used to compute mean values before opto onset
%     t_start_after_azi = 0.5*sr; %used to sort change in azimuth
    t_azi_end_loc = 0.25*sr; %250ms after opto to check if the eye moved
    t_azi_end_dur = 0.25*sr;
    t_sort_speed= 0.25*sr;
    t_sort=0.5*sr; 

    %pupil area vals analysis
    t_start_before_pa  = 0.25*sr;
    t_end_after_pa = 0.25*sr;

    %loco vals analysis
    t_start_before_l = 0.25*sr; %used to compute mean values before opto onset
    t_start_after_l = 0.5*sr; %sort by the mean value in this interval

elseif flash_dur_sec==0.2
    %% values for the analysis of 0.2sec flash
    flash_dur_sec=0.2;
    flash_dur=flash_dur_sec*sr;
    tbefore=flash_dur;
    tafter=2*flash_dur;

    %azimuth vals analysis
    t_start_before_azi = tbefore/2; %used to compute mean values before opto onset
    t_start_after_azi = flash_dur;

    %pupil vals analysis
    t_start_before_p = tbefore/2;

    %loco vals analysis
    t_start_before_l = tbefore; %used to compute mean values before opto onset
    t_start_after_l = flash_dur;
else
    warning('set analysis params');
    nf=0;
end

%% uncomment section below to process new recordings, otherwise will load below
% %% loop thru all the recordings
% % init vars
% fw_all=[];
% sw_all=[];
% rot_all=[];
% r_all=[];
% az_all=[];
% el_all=[];
% nr=0;
% rec_id=[];
% 
% for fi=1:nf
%     expfolder=folderlist{fi};    
%     folder_exp =  fullfile(repo_folder,expfolder);
%     save_file=fullfile(folder_exp ,'data.mat');
%     load(save_file); %should have data sturcture with trigger data now
% 
%     if data.stimparams.flash_t~=flash_dur_sec, continue;
%     else, disp(expfolder);
%     end
% 
% 
%     sr=data.trigger_data.sampling_rate;
%     repdur=data.stimparams.bout_len;
%     flash_dur_i = data.stimparams.flash_t;
%     start_opto_all= find(diff(data.trigger_data.opto_cont)==1)+1;    
%     nopto = length(start_opto_all);
%     start_opto=[];    
% 
%     for oi=1:nopto
%         %visual flash before opto
%         if any(data.trigger_data.red_signal(start_opto_all(oi)-tbefore:start_opto_all(oi)+flash_dur-1))
%             continue;
%         end
%         %other opto pulse too close before
%         if any(data.trigger_data.opto_cont(start_opto_all(oi)-tbefore:start_opto_all(oi)-1))
%             continue;
%         end
%         start_opto=[start_opto,start_opto_all(oi)];
%     end
%     nsi=length(start_opto);
% 
%     %% collect behav data
% 
%     %get continuous loco vals
%     loco_fw=interp1(data.trigger_data.loco.ball_time_s*sr,data.trigger_data.loco.forward_cm_s,1:length(data.trigger_data.opto_cont),"linear","extrap");
%     loco_sw=interp1(data.trigger_data.loco.ball_time_s*sr,data.trigger_data.loco.sideways_cm_s,1:length(data.trigger_data.opto_cont),"linear","extrap");
%     loco_rot=interp1(data.trigger_data.loco.ball_time_s*sr,data.trigger_data.loco.rotation_deg_s,1:length(data.trigger_data.opto_cont),"linear","extrap");
% 
%     tdur=tbefore+tafter+flash_dur;
%     tafterall=tafter+flash_dur;
%     minvals=min(length(data.trigger_data.cam_trigger), length(data.dlc.pupil_area));
% 
%     %remove outliers from eye values
%     r_clean=data.dlc.pupil_area((1:minvals));
%     idnan=isoutlier(r_clean,'movmedian',30);
%     r_clean(idnan)=NaN;
% 
%     el_clean=data.dlc.pupil_el((1:minvals));
% %     idnan=isoutlier(el_clean);
% %     el_clean(idnan)=NaN;
% 
%     az_clean=data.dlc.pupil_az((1:minvals));
% %     idnan=isoutlier(az_clean,'movmedian',30);
% % %     idnan=isoutlier(az_clean);
% %     az_clean(idnan)=NaN;
% 
%     %get continuous eye vals
%     r=interp1(data.trigger_data.cam_trigger(~idnan), r_clean(~idnan), 1:length(data.trigger_data.opto_cont),"linear","extrap");
%     az=interp1(data.trigger_data.cam_trigger(1:minvals), az_clean, 1:length(data.trigger_data.opto_cont),"linear","extrap");
%     el=interp1(data.trigger_data.cam_trigger(1:minvals), el_clean, 1:length(data.trigger_data.opto_cont),"linear","extrap");
% 
%     fw_o=zeros(nsi,tdur);
%     sw_o=zeros(nsi,tdur);
%     rot_o=zeros(nsi,tdur);
% 
%     r_o=zeros(nsi,tdur);
%     az_o=zeros(nsi,tdur);
%     el_o=zeros(nsi,tdur);
% 
%     for oi=1:nsi
%         fw_o(oi,:)=loco_fw(start_opto(oi)-tbefore:start_opto(oi)+tafterall-1);
%         sw_o(oi,:)=loco_sw(start_opto(oi)-tbefore:start_opto(oi)+tafterall-1);
%         rot_o(oi,:)=loco_rot(start_opto(oi)-tbefore:start_opto(oi)+tafterall-1);
% 
% 
%         r_o(oi,:)=r(start_opto(oi)-tbefore:start_opto(oi)+tafterall-1);
%         az_o(oi,:)=az(start_opto(oi)-tbefore:start_opto(oi)+tafterall-1);
%         el_o(oi,:)=el(start_opto(oi)-tbefore:start_opto(oi)+tafterall-1);
%     end
% 
%     fw_all=cat(1,fw_all,fw_o);
%     sw_all=cat(1,sw_all,sw_o);
%     rot_all=cat(1,rot_all,rot_o);
% 
%     r_all=cat(1,r_all,r_o);
%     az_all=cat(1,az_all,az_o);
%     el_all=cat(1,el_all,el_o);
%     nrep=size(r_o,1);
%     rec_id=cat(1,rec_id,ones(nrep,1)*fi);
% 
%     nr=nr+1;
%     behav_data(nr).data=save_file;
%     behav_data(nr).eye.r_on_opto=r_o;
%     behav_data(nr).eye.az_on_opto=az_o;
%     behav_data(nr).eye.el_on_opto=el_o;
%     behav_data(nr).loco.fw_opto=fw_o;
%     behav_data(nr).loco.sw_on_opto=sw_o;
%     behav_data(nr).loco.rot_on_opto=rot_o;
%     behav_data(nr).opto_trig_samples=start_opto;
% end
% 
% % remove outliers - errors of DLC
% maz=mean(az_all(:)); stdaz=std(az_all(:));
% mel=mean(el_all(:)); stdel=std(el_all(:));
% mr=mean(r_all(:)); stdr=std(r_all(:));
% idoutlier = unique([find(any(az_all>maz+5*stdaz,2)); find(any(az_all<maz-5*stdaz,2)); ...
%     find(any(el_all>mel+5*stdel,2)); find(any(el_all<mel-5*stdel,2)); find(any(r_all>mr+5*stdr,2)); find(any(r_all<mr-5*stdr,2))]);
% az_all(idoutlier,:)=[];
% el_all(idoutlier,:)=[];
% r_all(idoutlier,:)=[];
% rec_id(idoutlier,:)=[];
% 
% rec_id_u=unique(rec_id);
% rec_list=folderlist(rec_id);
% animal_id = unique(get_animal_id_from_recordingname(rec_list));
% 
% data_file=fullfile(path_to_figures, [file_str,'.mat']);
% save(data_file,'behav_data','fw_all', 'sw_all',...
%     'rot_all', 'r_all','az_all', 'el_all', 'rec_id','rec_list', 'animal_id');

%% load data
data_file=fullfile(path_to_figures, [file_str,'.mat']);
load(data_file);

%% plot

% sort_by_speed=1;
% figname=fullfile(path_to_figures,['az_change_velthresh_',num2str(round(flash_dur_sec*1000)),'ms']);
% scatter_delta_azimuth_on_start(az_all, tbefore, t_start_before_azi, t_azi_end_dur, t_azi_end_loc, t_sort, flash_dur, sort_by_speed, figname);

%make two plots selecting trials by speed
vel_thresh=0.7;
velstr=[num2str(fix(vel_thresh)),num2str(fix(vel_thresh*10)-fix(vel_thresh)*10)];
figname=fullfile(path_to_figures,[file_str,'_az_change_velthresh_',velstr]);
azimuth_plots_trials(az_all, tbefore, t_start_before_azi, t_azi_end_dur, t_azi_end_loc, t_sort_speed, flash_dur, vel_thresh,  rec_id, folderlist, figname);

%% raster plot of azimuth changes
figname=fullfile(path_to_figures,[file_str,'_az_change_tflash']);
make_raster_plot_azimuth_traces(az_all, tbefore, t_start_before_azi, t_sort,flash_dur,sr, rec_id, folderlist, figname);

%% plot pupil diameter changes
figname=fullfile(path_to_figures,[file_str,'_pupil_diam_tflash']);
make_pupil_area_plot(r_all, tbefore, t_start_before_pa, t_end_after_pa, flash_dur,sr, rec_id, folderlist, figname);

%% locomotion raster plots
figname=fullfile(path_to_figures,[file_str,'_loco_tflash']);
make_loco_plots(fw_all, sw_all, tbefore, t_start_before_l, t_sort, flash_dur,sr, rec_id, folderlist, figname);


%% loco videos upon stimulation


% %% raster plot of azimuth velocity
% figname=fullfile(path_to_figures,['az_vel_tflash',num2str(round(flash_dur_sec*1000)),'ms']);
% % make_raster_plot_azimuth_velocity(az_all, el_all, tbefore, t_start_before_azi, flash_dur, sr, figname);
% make_raster_plot_azimuth_velocity(az_all, tbefore, t_start_before_azi, t_sort, flash_dur, sr, figname);

% %% scater start and end of azimuth
% figname=fullfile(path_to_figures,['az_beforeafter_',num2str(round(flash_dur_sec*1000)),'ms']);
% sort_by_speed=1;
% % scatter_extreme_positions_before_after(az_all, el_all, tbefore, 0.25*sr, t_start_after_azi,flash_dur,sr, 100, figname);
% % scatter_start_end_azimuth(az_all, tbefore, t_start_before_azi, t_azi_end_dur, t_azi_end_loc, t_sort, figname);
% scatter_start_end_azimuth(az_all, tbefore, t_start_before_azi, t_azi_end_dur, t_azi_end_loc, t_sort, flash_dur, sort_by_speed, figname);
% sort_by_speed=0;
% scatter_start_end_azimuth(az_all, tbefore, t_start_before_azi, t_azi_end_dur, t_azi_end_loc, t_sort, flash_dur, sort_by_speed, figname);



% %% videos of extreme azimuth changes
% path_to_videos = 'D:\paper_figs\eye_opto_change';
% make_eye_azimuth_change_videos(behav_data,10, tbefore, t_start_before_azi, t_start_after_azi, flash_dur, path_to_videos);






function make_raster_plot_azimuth_traces(az_all, tbefore, tbefore_start, t_sort, flash_dur, sr, rec_id, rec_list, figname)

%get the number of reps, recordings and animals
nrep=size(az_all,1);
rec_id_u=unique(rec_id);
rlist=rec_list(rec_id_u);
nrec=numel(rlist);
animal_id = unique(get_animal_id_from_recordingname(rlist));
nanimals=numel(animal_id);
substr_n=['nrep=',num2str(nrep),', nrec=',num2str(nrec),', nanimals=',num2str(nanimals)];


az_start=mean(az_all(:,max(1,tbefore-tbefore_start):tbefore),2,"omitnan");
az_all_change=az_all-az_start;
%mean response to opto within the first tafter ms
meanaz=mean(az_all_change(:,tbefore+1:tbefore+t_sort),2);
[meanaz_sorted,sort_idk]=sort(meanaz);
% [~,sort_idk]=sort(az_start);


fig= figure; t=tiledlayout(1,12,"TileSpacing","compact"); cm=bluered;
ntraces=length(az_start);
nexttile([1,4]); % az starting position
az_start_sorted=az_start(sort_idk);
plot(az_start_sorted,1:ntraces,'k', LineWidth=1.5);
ylim([1,ntraces]);
yticks([]);
xticks([min(az_start(:)), max(az_start(:))]);
xticklabels(round([min(az_start(:)), max(az_start(:))],1));
ylabel('starting az');
ax1=gca;
set(ax1,'YDir','reverse');

subbin=20; %subsample every 5ms for plotting
nexttile([1,8]); %az change raster plot
maxv=max(abs(az_all_change(:)));
az_all_change_sorted=az_all_change(sort_idk,:);
imagesc(az_all_change_sorted(:,1:subbin:end),[-maxv,maxv]);
hold on;
plot([tbefore,tbefore]/subbin,[0.5,size(az_all_change,1)+0.5],'--k');
hold on;
plot([tbefore+flash_dur,tbefore+flash_dur]/subbin,[0.5,size(az_all_change,1)+0.5],'--k');
ax2=gca;
colormap(ax2, cm);
xticks([1,tbefore, tbefore+flash_dur]/subbin);
xticklabels(round([1,tbefore, tbefore+flash_dur]/sr,1));
xlabel('secs');
colorbar();

linkaxes([ax1,ax2],'y');

title(t,{'change in azimuth during opto'; substr_n});

% mean_trial_change = mean(abs(az_all));
% figure, plot(mean_trial_change);
% xticks([1,tbefore, tbefore+flash_dur]);
% xticklabels(round([1,tbefore, tbefore+flash_dur]/sr,1));

% varaz = var(az_all);
% figure, plot(varaz);
% xticks([1,tbefore, tbefore+flash_dur]);
% xticklabels(round([1,tbefore, tbefore+flash_dur]/sr,1));


savefig(fig, [figname,'.fig']);
% saveas(fig, [figname,'.svg']);
exportgraphics(fig, [figname,'.pdf'], 'BackgroundColor','none', 'ContentType','vector');
save([figname,'_data.mat'],'az_start_sorted', 'az_all_change_sorted', 'meanaz_sorted',...
    'tbefore', 'flash_dur',  'nrep', 'nrec', 'nanimals', 'rec_id', 'rec_list');
end




function make_eye_azimuth_change_videos(behav_data,nvid, tbefore, tbefore_start, tafter, flash_dur, path_to_videos)
%% find extreme azimuth changes during the first tafter ms
allmaxaz=[];
for s=1:length(behav_data)    
    ntrigs=length(behav_data(s).opto_trig_samples);
    for ti=1:ntrigs    
        azvals=behav_data(s).eye.az_on_opto(ti,:);
        azchange=azvals(tbefore+1:tbefore+flash_dur)-mean(azvals(max(1,tbefore-tbefore_start):tbefore));
        maxchange=mean(azchange(1:tafter));
        allmaxaz=[allmaxaz,maxchange];
    end
end
allmaxaz=sort(allmaxaz);
%take extreme nvid values
allmaxaz=[allmaxaz(1:nvid),allmaxaz(end-nvid+1:end)];

%find 10 extreme azimuth starting postions
% allmaxaz=[];
% for s=1:length(behav_data)    
%     ntrigs=length(behav_data(s).opto_trig_samples);
%     for ti=1:ntrigs
%         azvals=behav_data(s).eye.az_on_opto(ti,:);
%         azstart=mean(azvals(max(1,tbefore-tbefore_start):tbefore));
%         allmaxaz=[allmaxaz,azstart];
%     end
% end
% allmaxaz=sort(allmaxaz);
% allmaxaz=[allmaxaz(1:nvid),allmaxaz(end-nvid+1:end)];


%for the extreme eye movements, save videos
for s=1:length(behav_data)    
    ntrigs=length(behav_data(s).opto_trig_samples);
    for ti=1:ntrigs        
        azvals=behav_data(s).eye.az_on_opto(ti,:);
        azchange=azvals(tbefore+1:tbefore+flash_dur)-mean(azvals(max(1,tbefore-tbefore_start):tbefore));
        azi=mean(azchange(1:tafter));        
%         azi=mean(azvals(max(1,tbefore-tbefore_start):tbefore)); %starting position
        if ismember(azi,allmaxaz)
            %collect movie frames
            opto_sample=behav_data(s).opto_trig_samples(ti);
            load(behav_data(s).data);
            frnum = get_movie_frame_on_opto(data,opto_sample);
            %open the video
            oldavi_full=data.motion_energy.crop_avi;
            [~,oldavi_name,~]=fileparts(oldavi_full);
            if  ismember(azi,allmaxaz(1:nvid))  
                newavi=fullfile(fullfile(path_to_videos,'nasal'),[oldavi_name,'_',num2str(frnum),'_az',num2str(round(abs(azi))),'.avi']);
            else                
                newavi=fullfile(fullfile(path_to_videos,'temporal'),[oldavi_name,'_',num2str(frnum),'_az',num2str(round(abs(azi))),'.avi']);
            end            
            f0=frnum-30;
            fN=frnum+30;
            crop_video(oldavi_full, newavi, f0,fN);
        end
    end
end

end




function azimuth_plots_trials(az_all, tbefore, tbefore_start, tafter_start, toffset, tsort, flash_dur, vel_threshold,  rec_id, rec_list, figname)
startflash = tbefore_start+1;
endflash=startflash + flash_dur;
sr=20000;
ntrials=size(az_all,1);

%bin velocity and sort
%bin az values
binw=1/30*20000;
bine_s=round([1:binw:size(az_all,2)-1, size(az_all,2)]);
nbin=length(bine_s)-1;

[~,bin0_e]=min(abs(bine_s-(tbefore+1)));
[~,binN_e]=min(abs(bine_s-(tbefore+flash_dur)));

az_binned=zeros(size(az_all,1),nbin);
for bi=1:nbin
    az_binned(:,bi)=mean(az_all(:,bine_s(bi):bine_s(bi+1)),2);
end
%speed of az position as central difference
az_vel_binned=zeros(size(az_binned));
az_vel_binned(:,2:end-1)=(az_binned(:,3:end)-az_binned(:,1:end-2))/2;
az_vel_binned(:,1)=az_binned(:,2)-az_binned(:,1);
az_vel_binned(:,end)=az_binned(:,end)-az_binned(:,end-1);
az_vel_binned= az_vel_binned; %convert to change of az per sec

%find max velocity during flash and the width of the peak velocity
width_peaks=NaN(ntrials,1);
maxloc=NaN(ntrials,1);
maxvel=NaN(ntrials,1);
for i=1:ntrials
    [pks,locs,w,p] = findpeaks(abs(az_vel_binned(i,:)));
    foundpeak=0;    
    while foundpeak==0 && ~isempty(pks)
        [maxpeak,maxloc_id]=max(pks);
        % check if the peak happened during the flash
        if locs(maxloc_id)<bin0_e || locs(maxloc_id)>binN_e
            pks(maxloc_id)=[]; locs(maxloc_id)=[]; w(maxloc_id)=[];
            continue;
        end
        %check if the peak started before the flash
        if locs(maxloc_id)-w(maxloc_id)*2<bin0_e
            pks(maxloc_id)=[]; locs(maxloc_id)=[]; w(maxloc_id)=[];
            continue;
        end
        %keep this peak
        width_peaks(i)=w(maxloc_id);
        maxloc(i)=locs(maxloc_id);
        maxvel(i)=az_vel_binned(i,maxloc(i));
        foundpeak=1;
    end
end

%select trials with max vel above the threshold
sel_trials=find(abs(maxvel)>=vel_threshold);
speed_sel=az_vel_binned(sel_trials,:);
maxvel_sel=maxvel(sel_trials,:);
[speed_sorted,sort_idk]=sort(maxvel_sel);
%plot individual trial and check how to compute start/end position - maybe 
% using width?
az_start=arrayfun(@(i) az_binned(sel_trials(i),maxloc(sel_trials(i))-round(width_peaks(sel_trials(i))/2)-1), 1:numel(sel_trials));
az_end = arrayfun(@(i) az_binned(sel_trials(i),maxloc(sel_trials(i))+round(width_peaks(sel_trials(i))/2)+1), 1:numel(sel_trials));
% az_start=az_binned(sel_trials,maxloc(sel_trials)-round(width_peaks(sel_trials)/2)-1);
% az_end = az_binned(sel_trials,maxloc(sel_trials)+round(width_peaks(sel_trials)/2)+1);

% az_start=mean(az_binned(sel_trials,maxloc(sel_trials)-round(width_peaks(sel_trials)/2)-2:maxloc(sel_trials)-round(width_peaks(sel_trials)/2)-1),2);
% az_end = mean(az_binned(sel_trials,maxloc(sel_trials)+round(width_peaks(sel_trials)/2)+1:maxloc(sel_trials)+round(width_peaks(sel_trials)/2)+2),2);
speed_sorted_show=speed_sel(sort_idk,:);
az_start_sorted=az_start(sort_idk);
az_end_sorted=az_end(sort_idk);
daz_sorted = az_end_sorted - az_start_sorted;

%get the number of reps, recordings and animals of the selected trials
rec_id=rec_id(sel_trials);
nrep=numel(sel_trials);
rec_id_u=unique(rec_id);
rlist=rec_list(rec_id_u);
nrec=numel(rlist);
animal_id = unique(get_animal_id_from_recordingname(rlist));
nanimals=numel(animal_id);
substr_n=['nrep=',num2str(nrep),', nrec=',num2str(nrec),', nanimals=',num2str(nanimals)];


%bins for frequency counts of eye positions
edges_az=linspace(min(az_start_sorted), max(az_start_sorted),20);
binc_az=(edges_az(2)-edges_az(1))/2+edges_az(1:end-1);
edges_az_all=linspace(min(min(az_start_sorted),min(az_end_sorted)), max(max(az_start_sorted),max(az_end_sorted)),20);
binc_az_all=(edges_az_all(2)-edges_az_all(1))/2+edges_az_all(1:end-1);
edges_daz=linspace(min(daz_sorted), max(daz_sorted),20);
binc_daz=(edges_daz(2)-edges_daz(1))/2+edges_daz(1:end-1);

c_col=[143,46,119;53,143,141]/255;
c_col=cat(2,c_col,[0.3;0.3]);

st_az_1 = az_start_sorted(speed_sorted<=-vel_threshold);
st_az_2 = az_start_sorted(speed_sorted>=vel_threshold);
en_az_1 = az_end_sorted(speed_sorted<=-vel_threshold);
en_az_2 = az_end_sorted(speed_sorted>=vel_threshold);
daz_sorted_1= daz_sorted(speed_sorted<=-vel_threshold);
daz_sorted_2= daz_sorted(speed_sorted>=vel_threshold);

counts_az_s1 = histcounts(st_az_1,edges_az); %# get counts and bin locations
counts_az_s2 = histcounts(st_az_2,edges_az);
counts_daz_1 = histcounts(daz_sorted_1,edges_daz); %# get counts and bin locations
counts_daz_2 = histcounts(daz_sorted_2,edges_daz);

counts_az_s1_all = histcounts(st_az_1,edges_az_all); %# get counts and bin locations
counts_az_s2_all = histcounts(st_az_2,edges_az_all);
counts_en_s1_all = histcounts(en_az_1,edges_az_all); %# get counts and bin locations
counts_en_s2_all = histcounts(en_az_2,edges_az_all);

%% scatter plot of starting azimuth and change of az
f1= figure; t=tiledlayout(3,3);
nexttile([2,1]);

patch('XData', [0,counts_daz_1,0],'YData',[binc_daz(1), binc_daz, binc_daz(end)],...
    'FaceColor',c_col(1,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(1,1:3), 'EdgeAlpha', 0.5);
hold on;
patch('XData', [0,counts_daz_2,0],'YData',[binc_daz(1), binc_daz, binc_daz(end)],...
    'FaceColor',c_col(2,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(2,1:3), 'EdgeAlpha', 0.5);
ylabel('changed of  az');
ylim([binc_daz(1), binc_daz(end)]);
ax1=gca;

nexttile([2,2]);
scatter(st_az_1, daz_sorted_1, 15,'filled', 'MarkerFaceColor',c_col(1,1:3), 'MarkerFaceAlpha', 0.5);
hold on;
scatter(st_az_2, daz_sorted_2, 15,'filled', 'MarkerFaceColor',c_col(2,1:3), 'MarkerFaceAlpha', 0.5);
ylim([edges_daz(1), edges_daz(end)]);
xlim([edges_az(1), edges_az(end)]);
ax2=gca;

nexttile([1,1]);
axis off;

nexttile([1,2]);
patch('XData', [binc_az(1), binc_az, binc_az(end)],'YData',[0,counts_az_s1,0],...
    'FaceColor',c_col(1,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(1,1:3), 'EdgeAlpha', 0.5);
hold on;
patch('XData', [binc_az(1), binc_az, binc_az(end)],'YData',[0,counts_az_s2,0],...
    'FaceColor',c_col(2,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(2,1:3), 'EdgeAlpha', 0.5);
xlim([edges_az(1), edges_az(end)]);
ax3=gca;
xlabel('start az pos');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax3],'x');

title(t, {'scatter of change of azimuth vs starting position';['vel thresh = ',num2str(round(vel_threshold,1))]; substr_n});
figname1=[figname,'_1'];
savefig(f1,[figname1,'.fig']);
% saveas(f1,[figname1,'.svg']);
exportgraphics(f1, [figname1,'.pdf'], 'BackgroundColor','none', 'ContentType','vector');
saveas(f1,[figname1,'.png']);

%% scatter plot of start/end azimuth position
f1= figure; t=tiledlayout(3,3);
nexttile([2,1]);

patch('XData', [0,counts_en_s1_all,0],'YData',[binc_az_all(1), binc_az_all, binc_az_all(end)],...
    'FaceColor',c_col(1,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(1,1:3), 'EdgeAlpha', 0.5);
hold on;
patch('XData', [0,counts_en_s2_all,0],'YData',[binc_az_all(1), binc_az_all, binc_az_all(end)],...
    'FaceColor',c_col(2,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(2,1:3), 'EdgeAlpha', 0.5);
ylabel('end az');
ylim([edges_az_all(1), edges_az_all(end)]);
ax1=gca;

nexttile([2,2]);
scatter(st_az_1, en_az_1, 15,'filled', 'MarkerFaceColor',c_col(1,1:3), 'MarkerFaceAlpha', 0.5);
hold on;
scatter(st_az_2, en_az_2, 15,'filled', 'MarkerFaceColor',c_col(2,1:3), 'MarkerFaceAlpha', 0.5);
ylim([edges_az_all(1), edges_az_all(end)]);
xlim([edges_az_all(1), edges_az_all(end)]);
ax2=gca;

nexttile([1,1]);
axis off;

nexttile([1,2]);
patch('XData', [binc_az_all(1), binc_az_all, binc_az_all(end)],'YData',[0,counts_az_s1_all,0],...
    'FaceColor',c_col(1,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(1,1:3), 'EdgeAlpha', 0.5);
hold on;
patch('XData', [binc_az_all(1), binc_az_all, binc_az_all(end)],'YData',[0,counts_az_s2_all,0],...
    'FaceColor',c_col(2,1:3), 'FaceAlpha', 0.5,'EdgeColor',c_col(2,1:3), 'EdgeAlpha', 0.5);
xlim([edges_az_all(1), edges_az_all(end)]);
ax3=gca;
xlabel('start az pos');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax3],'x');

title(t, {'scatter of end vs starting azimuth position';['vel thresh = ',num2str(round(vel_threshold,1))]; substr_n});
figname1=[figname,'_3'];
savefig(f1,[figname1,'.fig']);
% saveas(f1,[figname1,'.svg']);
saveas(f1,[figname1,'.png']);
exportgraphics(f1, [figname1,'.pdf'], 'BackgroundColor','none', 'ContentType','vector');


%figure of velocity traces
fig= figure; t=tiledlayout(5,6,"TileSpacing","compact"); cm=bluered;


nexttile([3,4]);
maxv=max(abs(speed_sorted_show(:)));
imagesc(speed_sorted_show,[-maxv,maxv]);
hold on;
plot([bin0_e,bin0_e],[0.5,size(speed_sorted_show,1)+0.5],'--k');
hold on;
plot([binN_e,binN_e],[0.5,size(speed_sorted_show,1)+0.5],'--k');
ax1=gca;
colormap(ax1, cm);
xticks([1,bin0_e, binN_e]);
xticklabels(round([-tbefore, 0, flash_dur]/sr,2));
xlabel('secs');
colorbar();

nexttile([3,2]);
axis off;

nexttile([2,4]);
vel_mean=mean(abs(speed_sorted_show));
ntr=size(speed_sorted_show,1);
vel_stderr=std(abs(speed_sorted_show))/(ntr^0.5);
stdshade_mean_std_linecolor(vel_mean,vel_stderr,0.4,[0.5,0.5,0.5],[0,0,0]);
minvel=min(vel_mean-vel_stderr);
maxvel=max(vel_mean+vel_stderr);
hold on;
plot([bin0_e,bin0_e],[minvel,maxvel],'--k');
hold on;
plot([binN_e,binN_e],[minvel,maxvel],'--k');
ylim([minvel,maxvel]);
ax2=gca;
xticks([1,bin0_e, binN_e]);
xticklabels(round([-tbefore, 0, flash_dur]/sr,2));
linkaxes([ax1,ax2],'x');
xlim([1,size(speed_sorted_show,2)]);

nexttile([2,2]);
n_trials_selected=ntr;
n_trials_total=size(az_vel_binned,1);
n_low_vel=n_trials_total-n_trials_selected;
b= bar(1,[n_low_vel;n_trials_selected],'stacked');
b(1).FaceColor=[1,1,1];
b(2).FaceColor=[0.5,0.5,0.5];
ylim([0, n_trials_total+0.1*n_trials_total])
xticks([]);
yticks(n_low_vel);
yticklabels(round(n_low_vel/n_trials_total,2));
xlabel('% of sel trials');

title(t,{'pupil azimuth velocity';['vel thresh = ',num2str(round(vel_threshold,1))]; substr_n});
figname1=[figname,'_2'];
savefig(fig, [figname1,'.fig']);
saveas(fig, [figname1,'.png']);
% saveas(fig, [figname1,'.svg']);
exportgraphics(fig, [figname1,'.pdf'], 'BackgroundColor','none', 'ContentType','vector');

save([figname,'_data.mat'],...
    'st_az_1', 'st_az_2', 'daz_sorted_1', 'daz_sorted_2', ...    
    'counts_az_s1', 'counts_az_s2', 'counts_daz_1', 'counts_daz_2', ...
    'binc_az', 'binc_daz', 'tbefore', 'tbefore_start', 'tafter_start', 'toffset', 'tsort',...
    'binc_az_all','edges_az_all', 'counts_en_s1_all','counts_en_s2_all','counts_az_s1_all','counts_az_s2_all',...
    'st_az_1','en_az_1','st_az_2','en_az_2',...
    'speed_sorted_show', 'vel_mean','vel_stderr','bin0_e','binN_e','n_low_vel','n_trials_selected',...
     'width_peaks', 'maxloc', 'maxvel','az_vel_binned','az_binned',...
     'nrep', 'nrec', 'nanimals', 'rec_id', 'rec_list');
end


function make_pupil_area_plot(r_all, tbefore, tbefore_start, tafter_start, flash_dur, sr, rec_id, rec_list, figname)
% get the number of recordings and animals
nrep=size(r_all,1);
rec_id_u=unique(rec_id);
rlist=rec_list(rec_id_u);
nrec=numel(rlist);
animal_id = unique(get_animal_id_from_recordingname(rlist));
nanimals=numel(animal_id);
substr_n=['nrep=',num2str(nrep),', nrec=',num2str(nrec),', nanimals=',num2str(nanimals)];

%convert area_to_diameter
use_diam=1;
pix_mm=15; %measured the eye in the video
if use_diam
    d_all=((4*r_all/pi).^0.5)/pix_mm; %now in mm 
else
    d_all=r_all/(pix_mm*pix_mm); %now in mm^2 
end

d_before=mean(d_all(:,max(1,tbefore-tbefore_start):tbefore),2,"omitnan");
d_after =mean(d_all(:,tbefore+flash_dur:tbefore+flash_dur+tafter_start),2,"omitnan");
dmean=mean(d_all,1);
dstd=std(d_all,[],1);

%% test for the difference between d_before,d_after
% the pairs of [d_before(i),d_after(i)] are from the same trial, hence
% matched, use signed rank test
[test_res.p,test_res.h,test_res.stats] = signrank(d_before,d_after,'alpha',0.01);
test_res.name='signrank';

%% figure
if ~use_diam
    figname=strrep(figname,'diam','area');
end
fig=figure(Position=[100,100,500,500]); 
scatter(d_before,d_after,25,'filled', 'MarkerFaceColor',[0,0,0], 'MarkerFaceAlpha', 0.2);
hold on;
maxd=max([d_before',d_after']);
mind=min([d_before',d_after']);
plot([mind,maxd],[mind,maxd],'--k');
xticks([mind,maxd]);
xticklabels(round([mind,maxd],1));
yticks([mind,maxd]);
yticklabels(round([mind,maxd],1));
axis equal;
xlim([mind,maxd]);
ylim([mind,maxd]);
if use_diam
    xlabel('pupil diameter before opto (mm)');
    ylabel('pupil diameter after opto (mm)');
else
    xlabel('pupil area before opto (mm2)');
    ylabel('pupil area after opto (mm2)');
end
title({['d change, signrank test h=',num2str(test_res.h),', p=',num2str(test_res.p)]; substr_n});
savefig(fig,[figname,'_1.fig']);
% saveas(fig,[figname,'_1.svg']);
saveas(fig,[figname,'_1.pdf']);

scale_y=0.13;
fig=figure(Position=[100,100,500,500]); 
xvals=[1:size(r_all,2)]/sr;
ntraces2=size(r_all,1)^0.5;
if isempty(scale_y)
    ymin=min(dmean-dstd/ntraces2);
    ymax=max(dmean+dstd/ntraces2);
else
    ymin=mean(dmean)-scale_y/2;
    ymax=mean(dmean)+scale_y/2;
    yc=[mean(dmean)-1.1*scale_y/2, mean(dmean)+1.1*scale_y/2];
end


patch('XData',[(tbefore+1)/sr,(tbefore+1)/sr,(tbefore+flash_dur)/sr,(tbefore+flash_dur)/sr],...
    'YData',[ymin,ymax,ymax,ymin],...
     'FaceColor',[189, 229, 246]/255, 'FaceAlpha', 0.5,'EdgeColor',[189, 229, 246]/255, 'EdgeAlpha', 0.5);
hold on;
% options.handle=fig;
% options.color_area = [0.5,0.5,0.5];
% options.color_line = [0,0,0];
% options.alpha = 0.3;
%  options.line_width =1;
%  options.x_axis=xvals;
% options.error='sem';
% plot_areaerrorbar(d_all, options);

%subsample for the plot - every 5ms
dmean_s= dmean(1:20:end);
std_err_s=dstd(1:20:end)/ntraces2;
xvals_s=xvals(1:20:end);
stdshade_mean_std_linecolor(dmean_s,std_err_s,0.3,[0.5,0.5,0.5],[0,0,0],xvals_s);

if use_diam
    ylabel('pupil diameter (mm)');
else
    ylabel('pupil area (mm2)');
end
xlabel('trial duration (secs)');
if ~isempty(scale_y)
    ylim(yc);
end
title(substr_n);
savefig(fig,[figname,'_2.fig']);
% saveas(fig,[figname,'_2.svg']);
saveas(fig,[figname,'_2.pdf']);

save([figname,'_data.mat'], 'd_before', 'd_after', 'dmean', 'dstd', 'tbefore','flash_dur','test_res','ntraces2', 'nrep', 'nrec', 'nanimals', 'rec_id', 'rec_list');
end

function make_loco_plots(fw_all,sw_all, tbefore, tbefore_start, t_sort, flash_dur,sr, rec_id, rec_list, figname)
% make two raster plots for fw and sideways velocities
% sort the trials by the forward velocity in the first 500ms after opto onset

%tbefore is the time before opto onset in the traces fw_all and sw_all
% tbefore_start: time before opto used to compute starting loco pos

%get the number of reps, recordings and animals
nrep=size(fw_all,1);
rec_id_u=unique(rec_id);
rlist=rec_list(rec_id_u);
nrec=numel(rlist);
animal_id = unique(get_animal_id_from_recordingname(rlist));
nanimals=numel(animal_id);
substr_n=['nrep=',num2str(nrep),', nrec=',num2str(nrec),', nanimals=',num2str(nanimals)];


fw_before=mean(fw_all(:,max(1,tbefore-tbefore_start):tbefore),2);
fw_change=fw_all-fw_before;

fw_change=mean(fw_change(:,tbefore+1:tbefore+t_sort),2);
[~, idk_sorted]=sort(fw_change);
fw_all_sorted=fw_all(idk_sorted,:);
sw_all_sorted=sw_all(idk_sorted,:);

fig=figure; t=tiledlayout(1,2);
mval_fw=max(abs(fw_all_sorted(:)));
mval_sw=max(abs(sw_all_sorted(:)));
cm=bluered; 

%subsample for the plot
subbin=20; %sabsample every 5ms
nexttile;
imagesc(fw_all_sorted(:,1:subbin:end), [-mval_fw, mval_fw]);
hold on;
plot([tbefore+1,tbefore+1]/subbin,[0.5,nrep+0.5],'--', Color=[0.3,0.3,0.3]);
hold on;
plot([tbefore+flash_dur,tbefore+flash_dur]/subbin,[0.5,nrep+0.5],'--', Color=[0.3,0.3,0.3]);
xticks([1,tbefore,tbefore+flash_dur]/subbin);
xticklabels(round([1 tbefore tbefore+flash_dur]/sr,1));
ax1=gca;
colormap(ax1,cm);
hcb1 = colorbar;
hcb1.Title.String = "cm/sec";
title('forward velocity');

nexttile;
imagesc(sw_all_sorted(:,1:subbin:end), [-mval_sw, mval_sw]);
hold on;
plot([tbefore+1,tbefore+1]/subbin,[0.5,nrep+0.5],'--', Color=[0.3,0.3,0.3]);
hold on;
plot([tbefore+flash_dur,tbefore+flash_dur]/subbin,[0.5,nrep+0.5],'--', Color=[0.3,0.3,0.3]);
xticks([1,tbefore,tbefore+flash_dur]/subbin);
xticklabels(round([1 tbefore tbefore+flash_dur]/sr,1));
ax2=gca;
colormap(ax2,cm);
hcb2 = colorbar;
hcb2.Title.String = "cm/sec";
title('side velocity');
title(t, substr_n);

savefig(fig,[figname,'.fig']);
% saveas(fig,[figname,'.svg']);
exportgraphics(fig, [figname,'.pdf'], 'BackgroundColor','none', 'ContentType','vector');
save([figname,'_data.mat'], 'fw_all_sorted', 'sw_all_sorted', 'tbefore', 'flash_dur','sr', 'nrep', 'nrec', 'nanimals', 'rec_id', 'rec_list');
end

