%% script collects pupil values during black flashes
% The stimulus is 1 sec long flashes of various intencities (usually white,
% gray and black) randomly interleaved. The script extracts black 
% (from white to black) flashes. Plots: average pupil value per animal, slope,
% average pupil value while locomoting and still

% O.Symonova
% 2024

addpath 'D:\MATLAB\npxls\';
repof='D:\data\tnt_data';

resfolder = 'D:\data\resub\data\linear_step_flash';

% process_old_files='D:\data\tnt_data\linsteps_res\TNT_Sham_linflashres.mat';
% load(process_old_files);

%colldect all '_data.mat' files with linear flash stimulus
stimname='stimuli_linear_intensity';
datafilelist = collect_processed_stimuli_files(repof,stimname, 1, 0, 1);

%file to save results
linflash_datafile=fullfile(resfolder,'linflash_data.mat');

%parameters for the analysis
params_da.speed_thresh=5; %speed threshold 5cm/sec, below animal is still
params_da.tbefore_s=0.5; %time before flash onset
params_da.tafter_s=1; % time after flash onset
params_da.binw=0.025; % bin width in secs
params_da.pupil_area_pix2mmsq_scale=0.004; %(3^2mm^2/45^2pix) measurments from videos

binw=params_da.binw;
nbin=(params_da.tafter_s+params_da.tbefore_s)/binw; %number of bins
bine=linspace(-params_da.tbefore_s,params_da.tafter_s,nbin+1); %bin edges
histcen=bine(1:end-1) + binw/2; %bin centers
bin_st=interp1(histcen, 1:nbin, 0); %bin of the flash start

params_da.bin_st=bin_st;
params_da.histcen=histcen;

%% loop thru all recordings, separate all flash types, collect pupil responses
k=1;
linstep_pa_data={};
for fi=1:length(datafilelist)
    clear data;
    datafile=datafilelist{fi};
    load(datafile);

    if isfield(data.exp_info, 'folder')
        ofi=data.exp_info.folder;
    else
        ofi=data.exp_info.fullpath;
    end
    aniid = get_animal_id_from_recordingname(ofi);
    if contains(datafile,'control',"IgnoreCase",true), continue; end;

    [~,exp_name,~]= fileparts(datafile);
    if contains(exp_name,'tnt','IgnoreCase',true) %if the animal is TeLC or sham
        istnt=true;
    else
        istnt=false;
    end   

    %skip if triggers and camera frames do not match: alignment issue
    if  abs(numel(data.trigger_data.cam_trigger)-numel(data.dlc.pupil_az))>4
        continue;
    end   

    %find id of the logfile with flashes
    logid = find(contains({data.exp_info.logfiles(:).name},stimname));

    for li=1:length(logid)
        if isfield(data.exp_info.logfiles(logid), 'date')
            currdate=data.exp_info.logfiles(logid).date;
        else
            currdate=data.exp_info.logfiles(logid).file_tstamp;
        end
        stiparamfile = fullfile(ofi, data.exp_info.logfiles(logid(li)).name);
        stimdata = load(stiparamfile); %load stim params

        sr=data.trigger_data.sampling_rate;
        %start and end of the stimulus within the recording
        stimstart=data.exp_info.logfiles(logid).start_s;
        stimend=data.exp_info.logfiles(logid).end_s;

        %old recordings: cam triggers are saved in samples, convert to secs
        if sr<30000
            camtrig=data.trigger_data.cam_trigger/sr;
        else
            camtrig=data.trigger_data.cam_trigger;
        end
        %get all pupil values during the stimulus
        fr1=find(camtrig>=stimstart,1,"first");
        frn=find(camtrig<=stimend,1,"last");
        pa_all=data.dlc.pupil_area(fr1:frn);
        camtrig_stim=camtrig(fr1:frn);

        %get forward locomotion values during the stimulus
        bc1=find(data.trigger_data.loco.ball_time_s>=stimstart,1,"first");
        bcn=find(data.trigger_data.loco.ball_time_s<=stimend,1,"last");
        fwvals=data.trigger_data.loco.forward_cm_s(bc1:bcn);
        ballcamstim=data.trigger_data.loco.ball_time_s(bc1:bcn);

        %find mean speed during each frame -> intepolate speeds into frame bins
        speedframes=interp1(ballcamstim,fwvals,camtrig_stim,"linear");
        %get pupil values during still and moving periods independent from
        %flash types
        pa_still=pa_all(speedframes<=params_da.speed_thresh);
        pa_move=pa_all(speedframes>params_da.speed_thresh);

        % get pupil area during each flash

        %get the red signal during the recording and when each flash starts
        st=stimstart*sr;
        en=stimend*sr;
        stmin=max(1,round(st-params_da.tbefore_s*sr));
        enmax=min(length(data.trigger_data.red_signal), round(en+params_da.tbefore_s*sr));
        rs=data.trigger_data.red_signal(stmin:enmax);
        %start of each new flash
        vis_on = get_vis_on_edge_no_other(rs,zeros(size(rs)), 0.25, 0.25, sr);
        vis_on = (vis_on+double(st-params_da.tbefore_s*sr)-1)/sr;

        %there are 6 types of flashes: 0->1; 1->0; 0.5->1; 1->0.5; 0->0.5; 0.5->0;
        %from a random sequence of pairs of flashes create one vector:  0->1->0.5->1->0.5... 
        flash_id=nan(1,numel(stimdata.present_order_all)*2);
        nreps=size(stimdata.present_order_all,1);
        nlevels=size(stimdata.present_order_all,2);
        for i=1:nreps
            for j=1:nlevels
                ij=(i-1)*nlevels*2+(j-1)*2;
                pair=stimdata.present_order_all(i,j);
                flash_id(ij+1)=stimdata.pairsI(pair,1);
                flash_id(ij+2)=stimdata.pairsI(pair,2);
            end
        end

        
        nflash=length(flash_id);
       
        %for each flash type find reps
        flash_pairs=stimdata.pairsI;
        for pri=1:size(flash_pairs,1)
            idp =[];
            fpair=flash_pairs(pri,:);
            id1=find(flash_id==fpair(1)); %find occurences of starting intensity
            id1=id1(id1<nflash);
            for j=1:length(id1)
                if flash_id(id1(j)+1)==fpair(2) %if the next flash matches flash intensity in the pair
                    idp = [idp, id1(j)];
                end
            end

            flash_st=vis_on(idp)+1; %start of the flash in secs         

            %get the area of the pupil around the flash
            pupil_area=zeros(length(flash_st),nbin);           
            for vi=1:length(flash_st)
                triggi=camtrig(camtrig>flash_st(vi)-params_da.tbefore_s & camtrig<flash_st(vi)+params_da.tafter_s);
                f1=find(camtrig>(flash_st(vi)-params_da.tbefore_s),1,"first");
                f2=find(camtrig<(flash_st(vi)+params_da.tafter_s),1,"last");
                if f1>length(data.dlc.pupil_area) || f2>length(data.dlc.pupil_area)
                    pupil_area=pupil_area(1:vi-1,:);
                    break;
                end
                triggi=triggi-flash_st(vi);
                pupil_area(vi,:)=interp1(triggi,data.dlc.pupil_area(f1:f2),histcen,'linear','extrap');
            end
            linstep_pa_data(k).flash(pri).flash_st=stimdata.all_colors(fpair(1),:);
            linstep_pa_data(k).flash(pri).flash_en=stimdata.all_colors(fpair(2),:);
            linstep_pa_data(k).flash(pri).pupil_area=pupil_area;
        end
        linstep_pa_data(k).datafile=datafile;
        linstep_pa_data(k).aniid=aniid;
        linstep_pa_data(k).date=currdate;
        linstep_pa_data(k).istnt=istnt;
        linstep_pa_data(k).flash_dur_static=stimdata.tstatic;
        linstep_pa_data(k).flash_dur_lin=stimdata.t;
        linstep_pa_data(k).flash_colors=stimdata.all_colors;
        linstep_pa_data(k).stimstart=data.exp_info.logfiles(logid).start_s;
        linstep_pa_data(k).stimsend=data.exp_info.logfiles(logid).end_s;
        linstep_pa_data(k).pupil_area_all=pa_all;
        linstep_pa_data(k).pupil_area_all_sec =camtrig_stim;
        linstep_pa_data(k).loco_fw_all=fwvals;
        linstep_pa_data(k).loco_fw_all_sec=ballcamstim;
        linstep_pa_data(k).pupil_area_still=pa_still;
        linstep_pa_data(k).pupil_area_move=pa_move;
        k=k+1;
    end
end
save(linflash_datafile,'linstep_pa_data', 'params_da', '-v7.3');


%% find black flashes and pupil responses

load(linflash_datafile);
pupil_pix2mm=params_da.pupil_area_pix2mmsq_scale;
datevecall=sort(datetime({linstep_pa_data(:).date}))';
aniall=unique(cell2mat({linstep_pa_data(:).aniid})); %id of all animals

binw=params_da.binw;
nbin=(params_da.tafter_s+params_da.tbefore_s)/binw;
bine=linspace(-params_da.tbefore_s,params_da.tafter_s,nbin+1);
histcen=bine(1:end-1) + binw/2;
bin_st=interp1(histcen, 1:nbin, 0); %bin of the flash start
 

%% for each animal collect mean responses to the black flash, median pupil values
%values for animals from the control group
ctrl_traces=[];
ctrl_resp=[];
ctrl_slope=[];
ctrl_median_pa_still=[];
ctrl_median_pa_move=[];
ctrl_median_pa=[];
ctrl_all_pa_still=[];
ctrl_all_pa_move=[];
ctrl_all_pa=[];
ctrl_anilist=[];
%values for animals from the TeLC group
tnt_traces=[];
tnt_resp=[];
tnt_slope=[];  
tnt_median_pa_still=[];
tnt_median_pa_move=[];
tnt_median_pa=[];
tnt_all_pa_still=[];
tnt_all_pa_move=[];
tnt_all_pa=[];
tnt_anilist=[];
        
for ai=1:length(aniall)
    idai=find([linstep_pa_data(:).aniid]==aniall(ai));
    %take the first recording if multiple are available    
    [datevec_ani,ids]=sort(datetime({linstep_pa_data(idai).date}));
    ani_idone=idai(ids(1));    
    
    %find white-black flash
    wb_flash=find(...
        ismember(reshape([linstep_pa_data(ani_idone).flash(:).flash_st],3,[])',[0,255,255],'rows') &...
        ismember(reshape([linstep_pa_data(ani_idone).flash(:).flash_en],3,[])',[0,0,0],'rows')); %

    yvals=linstep_pa_data(ani_idone).flash(wb_flash).pupil_area;
    yvals_bl=yvals-mean(yvals(:,1:floor(bin_st)),2); %substract mean before flash
    yfull_resp=mean(yvals_bl); %mean response per animal

    bin200ms=find(histcen>=0.2,1,"first");
    yvals_resp=yfull_resp(bin200ms:end);
    xvals=histcen(bin200ms:end);
   
    %fit the line, use slope later
    [c,S] = polyfit(xvals,yvals_resp,1);

    %save full traces, y values used for the fit, and fitted slope
    if linstep_pa_data(ani_idone).istnt==1
        tnt_anilist=[tnt_anilist, aniall(ai)];
        tnt_traces=[tnt_traces; yfull_resp];        
        tnt_resp=[tnt_resp,yvals_resp];
        tnt_slope=[tnt_slope,c(1)];       
        tnt_median_pa_still=[tnt_median_pa_still, median(linstep_pa_data(ani_idone).pupil_area_still,'omitnan')];
        tnt_median_pa_move=[tnt_median_pa_move, median(linstep_pa_data(ani_idone).pupil_area_move,'omitnan')];
        tnt_median_pa=[tnt_median_pa, median(linstep_pa_data(ani_idone).pupil_area_all,'omitnan')];
        tnt_all_pa_still=[tnt_all_pa_still, linstep_pa_data(ani_idone).pupil_area_still'];
        tnt_all_pa_move=[tnt_all_pa_move,linstep_pa_data(ani_idone).pupil_area_move'];
        tnt_all_pa=[tnt_all_pa,linstep_pa_data(ani_idone).pupil_area_all'];
    else
        ctrl_anilist=[ctrl_anilist, aniall(ai)];
        ctrl_traces=[ctrl_traces; yfull_resp];        
        ctrl_resp=[ctrl_resp,yvals_resp];
        ctrl_slope=[ctrl_slope,c(1)];
        ctrl_median_pa_still=[ctrl_median_pa_still, median(linstep_pa_data(ani_idone).pupil_area_still,'omitnan')];
        ctrl_median_pa_move=[ctrl_median_pa_move, median(linstep_pa_data(ani_idone).pupil_area_move,'omitnan')];
        ctrl_median_pa=[ctrl_median_pa, median(linstep_pa_data(ani_idone).pupil_area_all,'omitnan')];
        ctrl_all_pa_still=[ctrl_all_pa_still, linstep_pa_data(ani_idone).pupil_area_still'];
        ctrl_all_pa_move=[ctrl_all_pa_move,linstep_pa_data(ani_idone).pupil_area_move'];
        ctrl_all_pa=[ctrl_all_pa,linstep_pa_data(ani_idone).pupil_area_all'];
    end         
end
%test the difference between the slopes
[ctrl_tnt_pupil_area.test_slope.pvalue, ctrl_tnt_pupil_area.test_slope.h0]=ranksum(ctrl_slope,tnt_slope);
ctrl_tnt_pupil_area.test_slope.test_type='ranksum';
%test the difference between the the median in still and move conditions
[ctrl_tnt_pupil_area.test_median_still.pvalue, ctrl_tnt_pupil_area.test_median_still.h0]=ranksum(ctrl_median_pa_still,tnt_all_pa_still);
ctrl_tnt_pupil_area.test_median_still.test_type='ranksum';
[ctrl_tnt_pupil_area.test_median_move.pvalue, ctrl_tnt_pupil_area.test_median_move.h0]=ranksum(ctrl_median_pa_move,tnt_all_pa_move);
ctrl_tnt_pupil_area.test_median_move.test_type='ranksum';
[ctrl_tnt_pupil_area.test_median_all.pvalue, ctrl_tnt_pupil_area.test_median_all.h0]=ranksum(ctrl_median_pa,tnt_all_pa);
ctrl_tnt_pupil_area.test_median_all.test_type='ranksum';

%% compute histogram of pupil values for plotting later
mina=min(min(tnt_all_pa_still),min(ctrl_all_pa_still));
maxa=max(max(tnt_all_pa_still),max(ctrl_all_pa_still));
histe_pa=linspace(mina,maxa,50);
histc_pa_still =histe_pa(1:end-1) + diff(histe_pa(1:2))/2;
histcount_ctrl_still=histcounts(ctrl_all_pa_still,histe_pa, "Normalization","probability");
histcount_tnt_still=histcounts(tnt_all_pa_still,histe_pa, "Normalization","probability");
mina=min(min(tnt_all_pa_move),min(ctrl_all_pa_move));
maxa=max(max(tnt_all_pa_move),max(ctrl_all_pa_move));
histe_pa=linspace(mina,maxa,50);
histc_pa_move =histe_pa(1:end-1) + diff(histe_pa(1:2))/2;
histcount_ctrl_move=histcounts(ctrl_all_pa_move,histe_pa, "Normalization","probability");
histcount_tnt_move=histcounts(tnt_all_pa_move,histe_pa, "Normalization","probability");
mina=min(min(tnt_all_pa),min(ctrl_all_pa));
maxa=max(max(tnt_all_pa),max(ctrl_all_pa));
histe_pa=linspace(mina,maxa,50);
histc_pa_all =histe_pa(1:end-1) + diff(histe_pa(1:2))/2;
histcount_ctrl_allpa=histcounts(ctrl_all_pa,histe_pa, "Normalization","probability");
histcount_tnt_allpa=histcounts(tnt_all_pa,histe_pa, "Normalization","probability");

%% save all the data
matfilename=fullfile(resfolder,'fullfieldflashes_pupil_values_data.mat');
save(matfilename,'ctrl_tnt_pupil_area','histcen','pupil_pix2mm',...
    'ctrl_traces', 'ctrl_resp', 'ctrl_slope', 'ctrl_median_pa_still','ctrl_median_pa_move', 'ctrl_median_pa',...
    'ctrl_all_pa_still','ctrl_all_pa_move', 'ctrl_all_pa', 'ctrl_anilist', ...
    'tnt_traces', 'tnt_resp', 'tnt_slope',   'tnt_median_pa_still', 'tnt_median_pa_move', 'tnt_median_pa',...
    'tnt_all_pa_still', 'tnt_all_pa_move', 'tnt_all_pa', 'tnt_anilist', ...
    'histc_pa_still','histcount_ctrl_still','histcount_tnt_still', ...
    'histc_pa_move','histcount_ctrl_move','histcount_tnt_move',...
    'histc_pa_all','histcount_ctrl_allpa','histcount_tnt_allpa',...
    '-v7.3');

%% pdf file where all the figure will be stored
figname='fullfieldflashes_pupil_values.pdf';
figname=fullfile(resfolder,figname);

%% figure for the full pupil change traces
fig = figure;
 plot(histcen,tnt_traces*pupil_pix2mm,'r');
 hold on; plot(histcen,ctrl_traces*pupil_pix2mm,'k');

maxval=max(max(ctrl_traces(:)), max(tnt_traces(:)));
minval=min(min(ctrl_traces(:)), min(tnt_traces(:)));
hold on; plot([0, 0],[minval*pupil_pix2mm,maxval*pupil_pix2mm],'--k');
xtickvals=round([histcen(1),0,histcen(end)],1);
xticks(xtickvals);
xlabel('secs');
ylabel('pupila area (mm^2)');
title('mean pupil change per animal');
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);

%% figure for the slopes
fig = figure;
bar([1:3],[mean(ctrl_slope)*pupil_pix2mm,NaN, mean(tnt_slope)*pupil_pix2mm],'FaceColor','none');
hold on;
xvals=[ones(1,length(ctrl_slope)), 3*ones(1,length(tnt_slope))];
s=swarmchart(xvals,[ctrl_slope*pupil_pix2mm,tnt_slope*pupil_pix2mm],'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'ctrl','tnt'});
ylabel('delta pupila area per sec (mm^2/sec)');
title({'pupil gradient'; ...
    ['ranksumtest: pvalue=',num2str(ctrl_tnt_pupil_area.test_slope.pvalue)]});
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);
   
%% figure for distribution of pupil area during the flashes - STILL
fig = figure;
plot(histc_pa_still*pupil_pix2mm,histcount_ctrl_still,'k'); hold on;
plot(histc_pa_still*pupil_pix2mm,histcount_tnt_still,'r'); hold on;
xlabel('pupil area - mm^2');
ylabel('prob.');
title('pupil area distribution STILL');
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);

%% figure for median pupil area during the flashes - STILL
fig = figure;
bar([1:3],[mean(ctrl_median_pa_still)*pupil_pix2mm,NaN, mean(tnt_median_pa_still)*pupil_pix2mm],'FaceColor','none');
hold on;
xvals=[ones(1,length(ctrl_median_pa_still)), 3*ones(1,length(tnt_median_pa_still))];
s=swarmchart(xvals,[ctrl_median_pa_still*pupil_pix2mm,tnt_median_pa_still*pupil_pix2mm],'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'ctrl','tnt'});
ylabel('pupil area - mm^2');
title({'median pupil area STILL'; ...
    ['ranksum pvalue=',num2str(ctrl_tnt_pupil_area.test_median_still.pvalue)]});
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);

%% figure for distribution of pupil area during the flashes - MOVE
fig = figure;
plot(histc_pa_move*pupil_pix2mm,histcount_ctrl_move,'k'); hold on;
plot(histc_pa_move*pupil_pix2mm,histcount_tnt_move,'r'); hold on;
xlabel('pupil area - mm^2');
ylabel('prob.');
title('pupil area distribution MOVE');
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);

%% figure for median pupil area during the flashes - MOVE
fig = figure;
bar([1:3],[mean(ctrl_median_pa_move)*pupil_pix2mm,NaN, mean(tnt_median_pa_move)*pupil_pix2mm],'FaceColor','none');
hold on;
xvals=[ones(1,length(ctrl_median_pa_move)), 3*ones(1,length(tnt_median_pa_move))];
s=swarmchart(xvals,[ctrl_median_pa_move*pupil_pix2mm,tnt_median_pa_move*pupil_pix2mm],'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'ctrl','tnt'});
ylabel('pupil area - mm^2');
title({'median pupil area MOVE'; ...
    ['ranksum pvalue=',num2str(ctrl_tnt_pupil_area.test_median_move.pvalue)]});
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);

%% figure for distribution of pupil area during the flashes - all values
fig = figure;
plot(histc_pa_all*pupil_pix2mm,histcount_ctrl_allpa,'k'); hold on;
plot(histc_pa_all*pupil_pix2mm,histcount_tnt_allpa,'r'); hold on;
xlabel('pupil area - mm^2');
ylabel('prob.');
title('pupil area distribution ALL');
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);

%% figure for median pupil area during the flashes - all values
fig = figure;
bar([1:3],[mean(ctrl_median_pa)*pupil_pix2mm,NaN, mean(tnt_median_pa)*pupil_pix2mm],'FaceColor','none');
hold on;
xvals=[ones(1,length(ctrl_median_pa)), 3*ones(1,length(tnt_median_pa))];
s=swarmchart(xvals,[ctrl_median_pa*pupil_pix2mm,tnt_median_pa*pupil_pix2mm],'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'ctrl','tnt'});
ylabel('pupil area - mm^2');
title({'median pupil area ALL'; ...
    ['ranksum pvalue=',num2str(ctrl_tnt_pupil_area.test_median_all.pvalue)]});
exportgraphics(fig,figname, 'BackgroundColor','none','ContentType','vector', 'Append',true);


