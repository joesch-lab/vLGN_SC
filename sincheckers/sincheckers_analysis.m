%% script collects all sinusoidal_checker stimuli
%for each recording the following values are collected:
% ani_id: animal id
% date: date of the recording;
% datafile: the processed recording data filr: triggers+spikes
% istnt: if the animal is TeLC
% is_wt: if the animal is WT
% is_npx: if the recording is Neuropixels
% bincen_t: time bins for eye values for each rep
% eye_vals_az: pupil azimuth in format [reps x bincen_t]
% eye_vals_el: same as above but elevation
% eye_vals_a: same as above but pupil area
% start_reps: time of the start of each repetition, aka new sine period
% moveyn: boolean vector classifying locomotion into still or moving
% camtrig_s: time stamp of frames of behavior video in secs


repof='D:\data\tnt_data';
resfolder='D:\data\resub\data\sin_checkers';


%colldect all '_data.mat' files
%check if they have stimuli_sinusoidal_checker stimulus
datafilelist = collect_processed_stimuli_files(repof,'stimuli_sinusoidal_checker', 1, 0, 1);

%% loop thru recordings, collect pupil and locomotion values
az_bincen_t=[0:0.1:4]; %each rep is 4sec, bin at 10ms
sindata={}; %structure to store extracted values
maxrep=Inf; %can specify how many reps to keep
k=1;
for fi=1:length(datafilelist)
    clear data;
    is_wt=false;

    datafile=datafilelist{fi};
    load(datafile);
    %if camera triggers and frames do not match: skip; alighnment problem
    if  abs(numel(data.trigger_data.cam_trigger)-numel(data.dlc.pupil_az))>4
        continue;
    end

    if isfield(data.exp_info, 'folder') %get full path of the experiment
        ofi=data.exp_info.folder;
    else
        ofi=data.exp_info.fullpath;
    end

    %check whether animal is TeLC, sham or WT
    if contains(ofi,'tnt','IgnoreCase',true)
        istnt=true;
    else
        istnt=false;
        if ~contains(ofi,'sham','IgnoreCase',true)
            is_wt=true;
        end
    end
    
    %get animal id from the recording
    [~,recnanme,~]=fileparts(ofi);
    ani_id= get_animal_id_from_recordingname({recnanme});    
    %skip wt 
    if is_wt
        continue;
    end

    %find id of the logfile with sincheckers
    logid = find(contains({data.exp_info.logfiles(:).name},'stimuli_sinusoidal_checker'));
    if isfield(data.exp_info.logfiles(logid), 'date')
        currdate=data.exp_info.logfiles(logid).date;
    else
        currdate=data.exp_info.logfiles(logid).file_tstamp;
    end

    sr=data.trigger_data.sampling_rate;

    %get start and end of the stimulus in samples 
    if isfield(data.exp_info.logfiles(logid),'start_s')
        st=data.exp_info.logfiles(logid).start_s*sr;
        en=data.exp_info.logfiles(logid).end_s*sr;
    else
        st=1;
        en=length(data.trigger_data.red_signal);
    end

    %red syncing signal: red signal is on whenever the stimulus changes
    %motion direction, find occurances of red signal to detect each rep
    rf=find(data.trigger_data.red_signal(st:en));
    drf=diff(rf); %time difference between red frames

    %find big gaps (each rep is 4 secs and includes motion to the left, 
    % then to the right; red frame indicates change of direction 
    % hence the gap btw red frames should be ~2sec)
    % a gap >2sec indicates delayed start of the stimulus, crop it out   
    hugegap=2*sr;
    hugegaploc=find(drf>hugegap);
    if ~isempty(hugegaploc)
        data.trigger_data.red_signal(1:rf(hugegaploc+1))=0;
        rf=find(data.trigger_data.red_signal);
        drf=diff(rf);
    end

    repgap=1*sr; %~2sec
    gaprep=find(drf>repgap);
    dirst=[rf(1), rf(gaprep+1)]; %start of new direction
    repst=dirst(1:2:end); %each rep has 2 directions
    repdur=min(diff(repst)); %duration of a rep
    repen=repst + repdur;
    nrep=length(repst);

    repst=(repst+st)/sr; %start of each rep in secs

    minfrnum = min(numel(data.trigger_data.cam_trigger), numel(data.dlc.pupil_az))-1;
    %conver triggers in old recordings to secs
    if data.trigger_data.sampling_rate>20000
        is_npx=1;        
        camtrig_s=data.trigger_data.cam_trigger(1:minfrnum);
    else
        is_npx=0;        
        camtrig_s=data.trigger_data.cam_trigger(1:minfrnum)/data.trigger_data.sampling_rate;
    end

    eyevals_az=data.dlc.pupil_az(1:minfrnum)';
    eyevals_el=data.dlc.pupil_el(1:minfrnum)';
    eyevals_a=data.dlc.pupil_area(1:minfrnum)';

    %for each rep find eye values
    allrep=min(nrep,maxrep);
    eye_vals_az=nan(allrep,numel(az_bincen_t));
    eye_vals_el=nan(allrep,numel(az_bincen_t));
    eye_vals_a=nan(allrep,numel(az_bincen_t));
    for ri=1:allrep
        %time vector for the current repetition
        xqt=repst(ri)+az_bincen_t;
        %find the az values at that time
        eye_vals_az(ri,:) = interp1(camtrig_s, eyevals_az, xqt, 'linear', 'extrap');
        eye_vals_el(ri,:) = interp1(camtrig_s, eyevals_el, xqt, 'linear', 'extrap');
        eye_vals_a(ri,:) = interp1(camtrig_s, eyevals_a, xqt, 'linear', 'extrap');
    end

    %get locomotion information
    avi_file = dir(fullfile(ofi,'*.avi'));
    avi_file=fullfile(avi_file.folder, avi_file.name);
    %still or moving by k-means clustering behav video
    move_yn=motion_from_video(avi_file);

    %save all data
    sindata(k).ani_id=ani_id;
    sindata(k).date=currdate;
    sindata(k).datafile=datafile;
    sindata(k).istnt=istnt;
    sindata(k).is_wt=is_wt;
    sindata(k).is_npx=is_npx;
    sindata(k).bincen_t=az_bincen_t;
    sindata(k).eye_vals_az=eye_vals_az;
    sindata(k).eye_vals_el=eye_vals_el;
    sindata(k).eye_vals_a=eye_vals_a;
    sindata(k).start_reps=repst;
    sindata(k).moveyn = move_yn;
    sindata(k).camtrig_s=camtrig_s;
    k=k+1;
end


T = struct2table(sindata); % convert the struct array to a table
sortedT = sortrows(T, 'date'); % sort the table by 'date'
sindata = table2struct(sortedT);
save(fullfile(resfolder, 'eye_vals_per_rep_all.mat'),'sindata');




%  %% for each recording construct the distrubution of spike numberes
% %% in the 100ms window for still/moving periods sham/tnt
% load(fullfile('D:\data\resub\data\sin_checkers', 'eye_vals_per_rep_all.mat'));
% figsavename=fullfile(resfolder, ['spike_countdistr_norm.pdf']);
% sample_dur=0.1; %s
% nsamples=1000;
% selrepnum=100;
% % spikecount_data={};
% % k=1;
% k=length(spikecount_data)+1;
% dat_inlist={spikecount_data(:).datafile};
% for fi=1:length(sindata)
%     clear data;
% 
%     if numel(sindata(fi).start_reps)<selrepnum
%         continue;
%     end
% 
%     currdate=sindata(fi).date(1:11);
%     datfile=sindata(fi).datafile;
% 
%     if ismember(datfile,dat_inlist)
%         continue;
%     end
% 
%     load(datfile);
%     if ~isfield(data, 'spike_data') || isempty(data.spike_data)
%         continue;
%     end
% 
%     %get spikes from all units
%     pop_spikes=cat(1,data.spike_data(:).spike_timing);
%     pop_spikes=sort(pop_spikes);
% 
%     if data.trigger_data.cam_trigger(2)>100
%         camtrig_s=data.trigger_data.cam_trigger/data.trigger_data.sampling_rate;
%         pop_spikes=double(pop_spikes)/data.trigger_data.sampling_rate;
%     else
%         camtrig_s=data.trigger_data.cam_trigger;
%     end
% 
%     minfrnum = min(numel(data.trigger_data.cam_trigger), numel(data.dlc.pupil_az));
%     camtrig_s=camtrig_s(1:minfrnum);
%     eyevals=data.dlc.pupil_area(1:minfrnum);
% 
% 
%     %get the on/off motion states
%     if isfield(data.exp_info, 'folder')
%         ofi=data.exp_info.folder;
%     else
%         ofi=data.exp_info.fullpath;
%     end
% 
%     avi_file = dir(fullfile(ofi,'*.avi'));
%     avi_file=fullfile(avi_file.folder, avi_file.name);
%     move_yn=motion_from_video(avi_file);
% 
%     stm=find(diff([0,move_yn'])==1);
%     enm=find(diff([move_yn',0])==-1);
%     t_st_motion=camtrig_s(stm+1);
%     t_en_motion=camtrig_s(enm);
% 
%     sts=find(diff([0,~move_yn'])==1);
%     ens=find(diff([~move_yn',0])==-1);
%     t_st_still=camtrig_s(sts+1);
%     t_en_still=camtrig_s(ens);
% 
%     %draw random samples
%     m_sample=sample_on_periods(t_st_motion, t_en_motion, sample_dur, nsamples);
%     s_sample=sample_on_periods(t_st_still, t_en_still, sample_dur, nsamples);
% 
%     m_bins=[m_sample; m_sample+sample_dur]';
%     m_bins=sortrows(m_bins,1);
%     s_bins=[s_sample; s_sample+sample_dur]';
%     s_bins=sortrows(s_bins,1);
% 
% 
%     %save spike counts for the random bins
%     spike_counts_ms=zeros(2,nsamples);
%     eyevals_insamples=zeros(2,nsamples);
%     for si=1:nsamples
%         %count spikes in moving bins
%         is=find(pop_spikes>=m_bins(si,1),1,"first");
%         ie=find(pop_spikes<=m_bins(si,2),1,"last");
%         if isempty(is) || isempty(ie) || ie-is<0
%             nc=0;
%         else
%             nc=ie-is+1;
%         end
%         spike_counts_ms(1,si)=nc;
%         %eye area in moving bins
%         camsel=camtrig_s>=m_bins(si,1) & camtrig_s<=m_bins(si,2);
%         eyevals_insamples(1,si)=mean(eyevals(camsel));
% 
%         %count spikes in still bins
%         is=find(pop_spikes>=s_bins(si,1),1,"first");
%         ie=find(pop_spikes<=s_bins(si,2),1,"last");
%         if isempty(is) || isempty(ie) || ie-is<0
%             nc=0;
%         else
%             nc=ie-is+1;
%         end
%         spike_counts_ms(2,si)=nc;
%         %eye area in still bins
%         camsel=camtrig_s>=s_bins(si,1) & camtrig_s<=s_bins(si,2);
%         eyevals_insamples(2,si)=mean(eyevals(camsel));
%     end
% 
%     %normalize population spike count by number of channels used in
%     %recording
%     if data.trigger_data.sampling_rate<21000 %npx has larger sampling rate
%         spike_counts_ms=spike_counts_ms/32;
%     else
%         spike_counts_ms=spike_counts_ms/384;
%     end
% 
%     %compute histograms of spike counts
%     min_spc=min(spike_counts_ms(:));
%     max_spc=max(spike_counts_ms(:));
%     spc_bine=linspace(min_spc,max_spc,51);
%     spc_binc=spc_bine(1:end-1)+diff(spc_bine)/2;
%     hc(1,:) = histcounts(spike_counts_ms(1,:),spc_bine,'Normalization','probability');
%     hc(2,:) = histcounts(spike_counts_ms(2,:),spc_bine,'Normalization','probability');
% 
%     %find the KL-divergence from still to moving distribution
%     distr_dist=0;
%     d_still=hc(2,:);
%     d_move=hc(1,:);
%     d_still=d_still+10e-4; %remove zeros
%     d_move=d_move+10e-4;
%     for bi=1:length(spc_binc)
%         distr_dist=distr_dist+d_move(bi)*log(d_move(bi)/d_still(bi));
%     end
%     % medians of ditstributions
%     median_still_move=median(spike_counts_ms,2);
%     % ks2 statistics
%     [ks_test.h,ks_test.p,ks_test.ks2stat] = kstest2(spike_counts_ms(2,:),spike_counts_ms(1,:));
%     ks_test.name='kstest2';
% 
% 
%     spikecount_data(k).ani_id = sindata(fi).ani_id;
%     spikecount_data(k).currdate = sindata(fi).date;
%     spikecount_data(k).datafile = sindata(fi).datafile;
%     spikecount_data(k).istnt = sindata(fi).istnt;
%     spikecount_data(k).is_wt = sindata(fi).is_wt;
%     spikecount_data(k).moveyn = move_yn;
%     spikecount_data(k).moving_bins = m_bins;
%     spikecount_data(k).still_bins = s_bins;
%     spikecount_data(k).spike_count_binc = spc_binc;
%     spikecount_data(k).spike_count_samples_moving_still=spike_counts_ms;
%     spikecount_data(k).spike_count_moving = hc(1,:);
%     spikecount_data(k).spike_count_still = hc(2,:);
%     spikecount_data(k).eyearea_moving_still=eyevals_insamples;
%     spikecount_data(k).distribution_medians_move_still=median_still_move;
%     spikecount_data(k).distribution_KLdist_still_move=distr_dist;
%     spikecount_data(k).distribution_ks2test=ks_test;
% 
% 
%     %make plots with the normalized spike count
%     titlestr0=['animal id: ',num2str(spikecount_data(k).ani_id), '; ',currdate];
%     if spikecount_data(k).istnt
%         titlestr0=[titlestr0, '; TNT'];
%     else
%         if spikecount_data(k).is_wt
%             titlestr0=[titlestr0, '; WT'];
%         else
%             titlestr0=[titlestr0, '; sham'];
%         end
%     end
%     move_dur=sum(move_yn)/length(move_yn);
%     titlestr1=['moving fraction: ', num2str(round(move_dur,3))];
%     titlestr2=['median diff: ', num2str(round(-diff(spikecount_data(k).distribution_medians_move_still),2)),...
%         '; KL-dist: ', num2str(round(spikecount_data(k).distribution_KLdist_still_move,2))];
%     titlestr3=[spikecount_data(k).distribution_ks2test.name, ...
%         ': pvalue=', num2str(round(spikecount_data(k).distribution_ks2test.p,4)),...
%         '; ks2stat= ', num2str(round(spikecount_data(k).distribution_ks2test.ks2stat,4))];
%     titlestr={titlestr0; titlestr1;titlestr2;titlestr3};
%     plot_spikecount_distr(spike_counts_ms(1,:), spike_counts_ms(2,:), 1, sample_dur, titlestr,figsavename);
% 
% 
%     k=k+1;
% end
% 
% save(fullfile(resfolder, 'spike_count_moving_still_normed.mat'),'spikecount_data');
% 


