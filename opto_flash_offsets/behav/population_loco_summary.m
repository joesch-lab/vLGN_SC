%% the script loads behav responses to opto flashes 
% from different groups and plots means

optodir='D:\data\resub\data\opto_locobehav\fig4';
ctrldir='D:\data\resub\data\opto_locobehav\sup_fig4';
ctrlSCinhdir='D:\data\resub\data\opto_locobehav\SC_inh_figs';

figsummary='D:\data\resub\data\opto_locobehav\opto_ctrl_summary.pdf';
figdata='D:\data\resub\data\opto_locobehav\opto_ctrl_summary_data.mat';

%% turning velocity difference
file_patt='optoloco*_loco_*data.mat';

datafile=dir(fullfile(optodir,file_patt));
datafile=fullfile(datafile.folder,datafile.name);
dataopto=load(datafile);

ctrlfile=dir(fullfile(ctrldir,file_patt));
ctrlfile=fullfile(ctrlfile.folder,ctrlfile.name);
datactrl=load(ctrlfile);

ctrlSCinhfile=dir(fullfile(ctrlSCinhdir,file_patt));
ctrlSCinhfile=fullfile(ctrlSCinhfile.folder,ctrlSCinhfile.name);
dataSCinhctrl=load(ctrlSCinhfile);

%time vector
binw=1/20000;
bine_s=-1:binw:2-binw/2;
bin_cen=bine_s(1:end-1)+diff(bine_s(1:2))/2;
nbin=length(bin_cen);
%find the bins of start opto and 250ms after the start
bin0=find(bine_s>=0,1,'first');
bin250ms=find(bine_s<=0.25,1,'last');
bin250msbefore=find(bine_s>=-0.25,1,'first');
%for each trial find the max velocity during [0,0.25]sec after opto
mean_before=mean(dataopto.sw_all_sorted(:,bin250msbefore:bin0),2);
maxval=max(dataopto.sw_all_sorted(:,bin0:bin250ms),[],2);
minval=min(dataopto.sw_all_sorted(:,bin0:bin250ms),[],2);
opto_maxswvel=choose_largest_extremum(minval, maxval)-mean_before;

mean_before=mean(datactrl.sw_all_sorted(:,bin250msbefore:bin0),2);
maxval=max(datactrl.sw_all_sorted(:,bin0:bin250ms),[],2);
minval=min(datactrl.sw_all_sorted(:,bin0:bin250ms),[],2);
ctrl_maxswvel=choose_largest_extremum(minval, maxval)-mean_before;

mean_before=mean(dataSCinhctrl.sw_all_sorted(:,bin250msbefore:bin0),2);
maxval=max(dataSCinhctrl.sw_all_sorted(:,bin0:bin250ms),[],2);
minval=min(dataSCinhctrl.sw_all_sorted(:,bin0:bin250ms),[],2);
ctrlSCinh_maxswvel=choose_largest_extremum(minval, maxval)-mean_before;

%% make figure of means and individual trials
mean_opto = mean(opto_maxswvel);
mean_ctrl = mean(ctrl_maxswvel);
mean_ctrlSCinh = mean(ctrlSCinh_maxswvel);
stderr_opto=std(opto_maxswvel)/(numel(opto_maxswvel)^0.5);
stderr_ctrl=std(ctrl_maxswvel)/(numel(ctrl_maxswvel)^0.5);
stderr_ctrlSCinh=std(ctrlSCinh_maxswvel)/(numel(ctrlSCinh_maxswvel)^0.5);
fig=figure;
bar([1,3,5], [mean_opto, mean_ctrl, mean_ctrlSCinh], 'FaceColor','none');
xticks([1,3,5]);
xticklabels({'opto', 'ctrl', 'inh SC'});
hold on;
plot([1,1],[mean_opto-stderr_opto,mean_opto+stderr_opto],'k',"Marker",".");
plot([3,3],[mean_ctrl-stderr_ctrl,mean_ctrl+stderr_ctrl],'k',"Marker",".");
plot([5,5],[mean_ctrlSCinh-stderr_ctrlSCinh,mean_ctrlSCinh+stderr_ctrlSCinh],'k',"Marker",".");
% x1=ones(size(opto_maxswvel))+rand(size(opto_maxswvel))*0.8-0.4;
% x2=ones(size(ctrl_maxswvel))*3+rand(size(ctrl_maxswvel))*0.8-0.4;
% x3=ones(size(ctrlSCinh_maxswvel))*5+rand(size(ctrlSCinh_maxswvel))*0.8-0.4;
% scatter([x1;x2; x3], [opto_maxswvel;ctrl_maxswvel; ctrlSCinh_maxswvel], '.k');
[~, pvalue_sw_velocity_opto_ctrl]=kstest2(opto_maxswvel,ctrl_maxswvel);
[~, pvalue_sw_velocity_opto_ctrlSCinh]=kstest2(opto_maxswvel,ctrlSCinh_maxswvel);
[~, pvalue_sw_velocity_ctrl_ctrlSCinh]=kstest2(ctrl_maxswvel,ctrlSCinh_maxswvel);
title({'max change in sideways velocity within 250ms';...
    ['kstest2: p(opto-ctrl)=',num2str(pvalue_sw_velocity_opto_ctrl),...
    '; p(opto-ctrlSCinh)=',num2str(pvalue_sw_velocity_opto_ctrlSCinh),'; p(ctrl-ctrlSCinh)=',num2str(pvalue_sw_velocity_ctrl_ctrlSCinh)];...
    ['opto: ntrials=', num2str(numel(opto_maxswvel)),'; nrec=',num2str(dataopto.nrec),'; nani=', num2str(dataopto.nanimals)]; ...
    ['ctrl: ntrials=', num2str(numel(ctrl_maxswvel)),'; nrec=',num2str(datactrl.nrec),'; nani=', num2str(datactrl.nanimals)];...
     ['ctrlSCinh: ntrials=', num2str(numel(ctrlSCinh_maxswvel)),'; nrec=',num2str(dataSCinhctrl.nrec),'; nani=', num2str(dataSCinhctrl.nanimals)]},...
     'FontSize',10);
exportgraphics(fig, figsummary, 'BackgroundColor','none', 'ContentType','vector', 'Append',true);
summary.sw_velocity.maxvals.opto=opto_maxswvel;
summary.sw_velocity.maxvals.ctrl=ctrl_maxswvel;
summary.sw_velocity.maxvals.ctrlSCinh=ctrlSCinh_maxswvel;
summary.sw_velocity.mean.opto=mean_opto;
summary.sw_velocity.mean.ctrl=mean_ctrl;
summary.sw_velocity.mean.ctrlSCinh=mean_ctrlSCinh;
summary.sw_velocity.stderr.opto=stderr_opto;
summary.sw_velocity.stderr.ctrl=stderr_ctrl;
summary.sw_velocity.stderr.ctrlSCinh=stderr_ctrlSCinh;
summary.sw_velocity.kstest2.pvalue_sw_velocity_opto_ctrl=pvalue_sw_velocity_opto_ctrl;
summary.sw_velocity.kstest2.pvalue_sw_velocity_opto_ctrlSCinh=pvalue_sw_velocity_opto_ctrlSCinh;
summary.sw_velocity.kstest2.pvalue_sw_velocity_ctrl_ctrlSCinh=pvalue_sw_velocity_ctrl_ctrlSCinh;

%% pupil diameter control
file_patt='optoloco*pupil_diam*data.mat';

datafile=dir(fullfile(optodir,file_patt));
datafile=fullfile(datafile.folder,datafile.name);
dataopto=load(datafile);

ctrlfile=dir(fullfile(ctrldir,file_patt));
ctrlfile=fullfile(ctrlfile.folder,ctrlfile.name);
datactrl=load(ctrlfile);

ctrlSCinhfile=dir(fullfile(ctrlSCinhdir,file_patt));
ctrlSCinhfile=fullfile(ctrlSCinhfile.folder,ctrlSCinhfile.name);
dataSCinhctrl=load(ctrlSCinhfile);

%delta of pupil diameter change
delta_diam_opto=dataopto.d_after-dataopto.d_before;
delta_diam_ctrl=datactrl.d_after-datactrl.d_before;
delta_diam_ctrlSCinh=dataSCinhctrl.d_after-dataSCinhctrl.d_before;


%make figure of means and individual trials
mean_opto = mean(delta_diam_opto);
mean_ctrl = mean(delta_diam_ctrl);
mean_ctrlSCinh = mean(delta_diam_ctrlSCinh);
stderr_opto=std(delta_diam_opto)/(numel(delta_diam_opto)^0.5);
stderr_ctrl=std(delta_diam_ctrl)/(numel(delta_diam_ctrl)^0.5);
stderr_ctrlSCinh=std(delta_diam_ctrlSCinh)/(numel(delta_diam_ctrlSCinh)^0.5);
fig=figure;
bar([1,3,5], [mean_opto, mean_ctrl, mean_ctrlSCinh], 'FaceColor','none');
xticks([1,3,5]);
xticklabels({'opto', 'ctrl', 'inh SC'});
hold on;
plot([1,1],[mean_opto-stderr_opto,mean_opto+stderr_opto],'k',"Marker",".");
plot([3,3],[mean_ctrl-stderr_ctrl,mean_ctrl+stderr_ctrl],'k',"Marker",".");
plot([5,5],[mean_ctrlSCinh-stderr_ctrlSCinh,mean_ctrlSCinh+stderr_ctrlSCinh],'k',"Marker",".");
[~, pvalue_deltadiam_opto_ctrl]=kstest2(delta_diam_opto,delta_diam_ctrl);
[~, pvalue_deltadiam_opto_ctrlSCinh]=kstest2(delta_diam_opto,delta_diam_ctrlSCinh);
[~, pvalue_deltadiam_ctrl_ctrlSCinh]=kstest2(delta_diam_ctrl,delta_diam_ctrlSCinh);
title({'delta pupil within 250ms';...
    ['kstest2: p(opto-ctrl)=',num2str(pvalue_deltadiam_opto_ctrl),...
    '; p(opto-ctrlSCinh)=',num2str(pvalue_deltadiam_opto_ctrlSCinh),'; p(ctrl-ctrlSCinh)=',num2str(pvalue_deltadiam_ctrl_ctrlSCinh)];...
    ['opto: ntrials=', num2str(numel(opto_maxswvel)),'; nrec=',num2str(dataopto.nrec),'; nani=', num2str(dataopto.nanimals)]; ...
    ['ctrl: ntrials=', num2str(numel(ctrl_maxswvel)),'; nrec=',num2str(datactrl.nrec),'; nani=', num2str(datactrl.nanimals)];...
     ['ctrlSCinh: ntrials=', num2str(numel(delta_diam_ctrlSCinh)),'; nrec=',num2str(dataSCinhctrl.nrec),'; nani=', num2str(dataSCinhctrl.nanimals)]},...
     'FontSize',10);
exportgraphics(fig, figsummary, 'BackgroundColor','none', 'ContentType','vector', 'Append',true);
summary.delta_diam_pupil.delta_vals.opto=delta_diam_opto;
summary.delta_diam_pupil.delta_vals.ctrl=delta_diam_ctrl;
summary.delta_diam_pupil.delta_vals.ctrlSCinh=delta_diam_ctrlSCinh;
summary.delta_diam_pupil.mean.opto=mean_opto;
summary.delta_diam_pupil.mean.ctrl=mean_ctrl;
summary.delta_diam_pupil.mean.ctrlSCinh=mean_ctrlSCinh;
summary.delta_diam_pupil.stderr.opto=stderr_opto;
summary.delta_diam_pupil.stderr.ctrl=stderr_ctrl;
summary.delta_diam_pupil.stderr.ctrlSCinh=stderr_ctrlSCinh;
summary.delta_diam_pupil.kstest2.pvalue_deltadiam_opto_ctrl=pvalue_deltadiam_opto_ctrl;
summary.delta_diam_pupil.kstest2.pvalue_deltadiam_opto_ctrlSCinh=pvalue_deltadiam_opto_ctrlSCinh;
summary.delta_diam_pupil.kstest2.pvalue_deltadiam_ctrl_ctrlSCinh=pvalue_deltadiam_ctrl_ctrlSCinh;

%% pupil velocity control
file_patt='optoloco*az_change_velthresh*.mat';

datafile=dir(fullfile(optodir,file_patt));
datafile=fullfile(datafile.folder,datafile.name);
dataopto=load(datafile);

ctrlfile=dir(fullfile(ctrldir,file_patt));
ctrlfile=fullfile(ctrlfile.folder,ctrlfile.name);
datactrl=load(ctrlfile);

ctrlSCinhfile=dir(fullfile(ctrlSCinhdir,file_patt));
ctrlSCinhfile=fullfile(ctrlSCinhfile.folder,ctrlSCinhfile.name);
dataSCinhctrl=load(ctrlSCinhfile);

%time vector
binw=1/30;
bine_s=-1:binw:2;
bin_cen=bine_s(1:end-1)+diff(bine_s(1:2))/2;
nbin=length(bin_cen);

%find the bins of start opto and 250ms after the start
bin0=find(bine_s>=0,1,'first');
bin250ms=find(bine_s<=0.25,1,'last');
%for each trial find the max velocity during [0,0.25]sec after opto
opto_maxpupvel=max(abs(dataopto.speed_sorted_show(:,bin0:bin250ms)),[],2);
ctrl_maxpupvel=max(abs(datactrl.speed_sorted_show(:,bin0:bin250ms)),[],2);
ctrlSCinh_maxpupvel=max(abs(dataSCinhctrl.speed_sorted_show(:,bin0:bin250ms)),[],2);

%make figure of means and individual trials
mean_opto = mean(opto_maxpupvel);
mean_ctrl = mean(ctrl_maxpupvel);
mean_ctrlSCinh = mean(ctrlSCinh_maxpupvel);
stderr_opto=std(opto_maxpupvel)/(numel(opto_maxpupvel)^0.5);
stderr_ctrl=std(ctrl_maxpupvel)/(numel(ctrl_maxpupvel)^0.5);
stderr_ctrlSCinh=std(ctrlSCinh_maxpupvel)/(numel(ctrlSCinh_maxpupvel)^0.5);
fig=figure;
bar([1,3,5], [mean_opto, mean_ctrl, mean_ctrlSCinh], 'FaceColor','none');
xticks([1,3,5]);
xticklabels({'opto', 'ctrl', 'inh SC'});
hold on;
plot([1,1],[mean_opto-stderr_opto,mean_opto+stderr_opto],'k',"Marker",".");
plot([3,3],[mean_ctrl-stderr_ctrl,mean_ctrl+stderr_ctrl],'k',"Marker",".");
plot([5,5],[mean_ctrlSCinh-stderr_ctrlSCinh,mean_ctrlSCinh+stderr_ctrlSCinh],'k',"Marker",".");
ylabel('deg/frame, 30fps');
[~, pvalue_pupvelocity_opto_ctrl]=kstest2(opto_maxpupvel,ctrl_maxpupvel);
[~, pvalue_pupvelocity_opto_ctrlSCinh]=kstest2(opto_maxpupvel,ctrlSCinh_maxpupvel);
[~, pvalue_pupvelocity_ctrl_ctrlSCinh]=kstest2(ctrl_maxpupvel,ctrlSCinh_maxpupvel);
title({'max pupil velocity within 250ms';...
    ['kstest2: p(opto-ctrl)=',num2str(pvalue_pupvelocity_opto_ctrl),...
    '; p(opto-ctrlSCinh)=',num2str(pvalue_pupvelocity_opto_ctrlSCinh),'; p(ctrl-ctrlSCinh)=',num2str(pvalue_pupvelocity_ctrl_ctrlSCinh)];...
    ['opto: ntrials=', num2str(numel(opto_maxpupvel)),'; nrec=',num2str(dataopto.nrec),'; nani=', num2str(dataopto.nanimals)]; ...
    ['ctrl: ntrials=', num2str(numel(ctrl_maxpupvel)),'; nrec=',num2str(datactrl.nrec),'; nani=', num2str(datactrl.nanimals)];...
     ['ctrlSCinh: ntrials=', num2str(numel(ctrlSCinh_maxpupvel)),'; nrec=',num2str(dataSCinhctrl.nrec),'; nani=', num2str(dataSCinhctrl.nanimals)]},...
     'FontSize',10);
exportgraphics(fig, figsummary, 'BackgroundColor','none', 'ContentType','vector', 'Append',true);
summary.pupil_maxvelocity.maxvals.opto=opto_maxpupvel;
summary.pupil_maxvelocity.maxvals.ctrl=ctrl_maxpupvel;
summary.pupil_maxvelocity.maxvals.ctrlSCinh=ctrlSCinh_maxpupvel;
summary.pupil_maxvelocity.mean.opto=mean_opto;
summary.pupil_maxvelocity.mean.ctrl=mean_ctrl;
summary.pupil_maxvelocity.mean.ctrlSCinh=mean_ctrlSCinh;
summary.pupil_maxvelocity.stderr.opto=stderr_opto;
summary.pupil_maxvelocity.stderr.ctrl=stderr_ctrl;
summary.pupil_maxvelocity.stderr.ctrlSCinh=stderr_ctrlSCinh;
summary.pupil_maxvelocity.kstest2.pvalue_pupvelocity_opto_ctrl=pvalue_pupvelocity_opto_ctrl;
summary.pupil_maxvelocity.kstest2.pvalue_pupvelocity_opto_ctrlSCinh=pvalue_pupvelocity_opto_ctrlSCinh;
summary.pupil_maxvelocity.kstest2.pvalue_pupvelocity_ctrl_ctrlSCinh=pvalue_pupvelocity_ctrl_ctrlSCinh;



%% pupil velocity control - cumulative eye motions
file_patt='optoloco*az_change_velthresh*.mat';

datafile=dir(fullfile(optodir,file_patt));
datafile=fullfile(datafile.folder,datafile.name);
dataopto=load(datafile);

ctrlfile=dir(fullfile(ctrldir,file_patt));
ctrlfile=fullfile(ctrlfile.folder,ctrlfile.name);
datactrl=load(ctrlfile);

ctrlSCinhfile=dir(fullfile(ctrlSCinhdir,file_patt));
ctrlSCinhfile=fullfile(ctrlSCinhfile.folder,ctrlSCinhfile.name);
dataSCinhctrl=load(ctrlSCinhfile);

%time vector
binw=1/30;
bine_s=-1:binw:2;
bin_cen=bine_s(1:end-1)+diff(bine_s(1:2))/2;
nbin=length(bin_cen);

%find the bins of start opto and 250ms after the start
bin0=find(bine_s>=0,1,'first');
bin250ms=find(bine_s<=0.25,1,'last');
%for each trial find the max velocity during [0,0.25]sec after opto
% opto_maxpupvel=sum(abs(dataopto.speed_sorted_show(:,bin0:bin250ms)),2);
% ctrl_maxpupvel=sum(abs(datactrl.speed_sorted_show(:,bin0:bin250ms)),2);
% ctrlSCinh_maxpupvel=sum(abs(dataSCinhctrl.speed_sorted_show(:,bin0:bin250ms)),2);
opto_maxpupvel=sum(dataopto.speed_sorted_show(:,bin0:bin250ms),2);
ctrl_maxpupvel=sum(datactrl.speed_sorted_show(:,bin0:bin250ms),2);
ctrlSCinh_maxpupvel=sum(dataSCinhctrl.speed_sorted_show(:,bin0:bin250ms),2);

%make figure of means and individual trials
mean_opto = mean(opto_maxpupvel);
mean_ctrl = mean(ctrl_maxpupvel);
mean_ctrlSCinh = mean(ctrlSCinh_maxpupvel);
stderr_opto=std(opto_maxpupvel)/(numel(opto_maxpupvel)^0.5);
stderr_ctrl=std(ctrl_maxpupvel)/(numel(ctrl_maxpupvel)^0.5);
stderr_ctrlSCinh=std(ctrlSCinh_maxpupvel)/(numel(ctrlSCinh_maxpupvel)^0.5);
fig=figure;
bar([1,3,5], [mean_opto, mean_ctrl, mean_ctrlSCinh], 'FaceColor','none');
xticks([1,3,5]);
xticklabels({'opto', 'ctrl', 'inh SC'});
hold on;
plot([1,1],[mean_opto-stderr_opto,mean_opto+stderr_opto],'k',"Marker",".");
plot([3,3],[mean_ctrl-stderr_ctrl,mean_ctrl+stderr_ctrl],'k',"Marker",".");
plot([5,5],[mean_ctrlSCinh-stderr_ctrlSCinh,mean_ctrlSCinh+stderr_ctrlSCinh],'k',"Marker",".");
ylabel('deg in 0.25sec');
[~, pvalue_pupshift_opto_ctrl]=kstest2(opto_maxpupvel,ctrl_maxpupvel);
[~, pvalue_pupshift_opto_ctrlSCinh]=kstest2(opto_maxpupvel,ctrlSCinh_maxpupvel);
[~, pvalue_pupshift_ctrl_ctrlSCinh]=kstest2(ctrl_maxpupvel,ctrlSCinh_maxpupvel);
title({'cumulative pupil shift within 250ms';...
    ['kstest2: p(opto-ctrl)=',num2str(pvalue_pupshift_opto_ctrl),...
    '; p(opto-ctrlSCinh)=',num2str(pvalue_pupshift_opto_ctrlSCinh),'; p(ctrl-ctrlSCinh)=',num2str(pvalue_pupshift_ctrl_ctrlSCinh)];...
    ['opto: ntrials=', num2str(numel(opto_maxpupvel)),'; nrec=',num2str(dataopto.nrec),'; nani=', num2str(dataopto.nanimals)]; ...
    ['ctrl: ntrials=', num2str(numel(ctrl_maxpupvel)),'; nrec=',num2str(datactrl.nrec),'; nani=', num2str(datactrl.nanimals)];...
     ['ctrlSCinh: ntrials=', num2str(numel(ctrlSCinh_maxpupvel)),'; nrec=',num2str(dataSCinhctrl.nrec),'; nani=', num2str(dataSCinhctrl.nanimals)]},...
     'FontSize',10);
fig_cumshift='D:\data\resub\data\opto_locobehav\opto_ctrl_pupil_shift250ms.pdf';
exportgraphics(fig, fig_cumshift, 'BackgroundColor','none', 'ContentType','vector', 'Append',true);

pupvel_opto=dataopto.speed_sorted_show;
pupvel_ctrl=datactrl.speed_sorted_show;
pupvel_ctrlSCinh=dataSCinhctrl.speed_sorted_show;
save('D:\data\resub\data\opto_locobehav\pupil_velocities.mat', 'bine_s', 'bin_cen','pupvel_opto', 'pupvel_ctrl', 'pupvel_ctrlSCinh', '-v7.3');

save(figdata,'summary', '-v7.3');



function extrvals= choose_largest_extremum(minvals, maxvals)
extrvals=zeros(size(minvals));
for i=1:length(minvals)
if abs(maxvals(i))>=abs(minvals(i))
    extrvals(i)=maxvals(i);    
else
    extrvals(i)=minvals(i);
end
end
end




