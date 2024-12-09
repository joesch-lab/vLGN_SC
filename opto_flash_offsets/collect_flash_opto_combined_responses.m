%this script collects responses (PSTH, spike count, time of the first spike)
%during visual flashes, opto stimulattion and their overlap

load('D:\data\resub\data\animal_ids'); % animal_ids
align_peaks=false;


ofresfolder = 'D:\data\resub\data';
list_process = 'D:\data\resub\data\with_laser_list.txt'; %laser_yn=1;
folder_list = get_list_from_file(list_process);
path_to_figures= 'D:\data\resub\figs';

pval_thresh=0.01;
arrange_by_dI=true;
all_responses=[];
all_params=[];
all_spikecount=[];
all_firstspike=[];
rec_id=[];
minD=0;
maxD=Inf;
splitD=500;
for di=1:length(folder_list)
    ofdata_file=fullfile(ofresfolder ,folder_list{di});
    disp([num2str(di), ' ', ofdata_file]);
    ofdata_file = fullfile(ofdata_file,'opto_flash_data.mat');
    [all_responses_i, max_fr_params, spike_count, first_spike, bine, sr] = get_firing_props_per_offset(ofdata_file, pval_thresh,arrange_by_dI);           
    all_responses=cat(1,all_responses,all_responses_i);
    all_params=cat(1,all_params,max_fr_params);
    all_spikecount=cat(1,all_spikecount,spike_count);
    all_firstspike=cat(1,all_firstspike,first_spike);
    rec_id=[rec_id, ones(1,size(all_responses_i,1))*di];   
end

disp('Done');
save(fullfile(path_to_figures, 'flash_opto_all_raw.mat'),'all_responses','all_params', 'all_spikecount','all_firstspike', 'bine', 'pval_thresh', 'rec_id', 'folder_list');

%% figure all responses for all neurons
figname = fullfile(path_to_figures,'flop_all_visopto_neurons');
make_raster_all_responses(all_responses, rec_id, folder_list, animal_ids, bine, figname);

%% figure for flash, flash+opto, opto responses for sSC
minD=0;
maxD=500;
figname = fullfile(path_to_figures,'flop_visopto_neurons_sSC');
make_raster_flon_flopon_opto(all_responses, rec_id, folder_list, animal_ids, bine, minD, maxD, 'sSC', figname);

%% figure for flash, flash+opto, opto responses for dSC
minD=500;
maxD=Inf;
figname = fullfile(path_to_figures,'flop_visopto_neurons_dSC');
make_raster_flon_flopon_opto(all_responses, rec_id, folder_list, animal_ids, bine, minD, maxD, 'dSC', figname);

% %% spike count plots for flash vs flash+opto
% make_spike_count_flash_opto_scatterplot(all_spikecount, rec_id, folder_list, animal_ids, 0,500,'sSC', path_to_figures);
% make_spike_count_flash_opto_scatterplot(all_spikecount, rec_id, folder_list, animal_ids, 500,Inf,'dSC', path_to_figures);
% % %make_spike_count_flash_opto_scatterplot(all_spikecount, rec_id, folder_list, animal_ids, 0,Inf,'all', path_to_figures);
% 
% % %% spike count plots for baseline vs baseline+opto
% make_spike_count_bl_blopto_scatterplot(all_spikecount, rec_id, folder_list, animal_ids, 0,500,'sSC', path_to_figures);
% make_spike_count_bl_blopto_scatterplot(all_spikecount, rec_id, folder_list, animal_ids, 500,Inf,'dSC', path_to_figures);
% % % make_spike_count_bl_blopto_scatterplot(all_spikecount, rec_id, folder_list, animal_ids, 0,Inf,'all', path_to_figures);
% 
% % time of the first spike for flash vs flash+opto
% make_first_spike_timing_flash_opto_scatterplot(all_firstspike, rec_id, folder_list, animal_ids, 0,500,'sSC', path_to_figures);
% make_first_spike_timing_flash_opto_scatterplot(all_firstspike, rec_id, folder_list, animal_ids, 500,Inf,'dSC', path_to_figures);
% % make_first_spike_timing_flash_opto_scatterplot(all_firstspike, rec_id, folder_list, animal_ids, 0,Inf,'all', path_to_figures);



function make_raster_flon_flopon_opto(all_responses, rec_id, rec_names, animal_ids, bine, minD, maxD, titlestr, figname)

idx=find(all_responses(:,2)>=minD & all_responses(:,2)<maxD);
rec_id=rec_id(idx);
responses_raw = all_responses(idx,3:end);
depth=all_responses(idx,2);
responses_raw = responses_raw(:,[1*60:1*60+60-1, 3*60:3*60+60-1, 4*60:4*60+60-1]); 
%normalize to max
max_n = max(responses_raw,[],2);
responses_n=responses_raw./max_n;
%sort time to maxFR
[~,t2max]=max(responses_n(:,1:60),[],2);
[~,idx]=sort(t2max);
responses_n=responses_n(idx,:);
max_n=max_n(idx);
depth=depth(idx);

nbin = length(bine)-1;
onset_i=find(bine>=0,1,"first");
flash_resp=responses_n(:,1:nbin);
opto_resp=responses_n(:,1+2*nbin:end);
flashopto_resp=responses_n(:,1+nbin:2*nbin);

f = figure; t=tiledlayout(5,3,"TileSpacing",'compact');
nn=size(responses_n,1);
t_ticks=1:10:nbin;
xticks_labels=bine(t_ticks);

nexttile([4,1]);
imagesc(flash_resp, [0,1]); title('vis off');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks([]);
xticks([]);
% xticks(t_ticks);
% xticklabels(xticks_labels);
ax1=gca;

nexttile([4,1]);
imagesc(opto_resp, [0,1]); title('opto');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks([]);
xticks([]);
% xticks(t_ticks);
% xticklabels(xticks_labels);
ax2=gca;

nexttile([4,1]);
imagesc(flashopto_resp, [0,1]); title('vis off + opto');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks([]);
xticks([]);
% xticks(t_ticks);
% xticklabels(xticks_labels);
colormap(flipud(gray));
ax3=gca;

%compute mean and stderr of normailized responses
mean_flash=mean(flash_resp,'omitnan');
std_flash=std(flash_resp,'omitnan');
mean_opto=mean(opto_resp,'omitnan');
std_opto=std(opto_resp,'omitnan');
mean_flashopto=mean(flashopto_resp,'omitnan');
std_flashopto=std(flashopto_resp,'omitnan');
maxval=max([max(mean_flash+std_flash), max(mean_opto+std_opto), max(mean_flashopto+std_flashopto)]);
minval=min([min(mean_flash-std_flash), min(mean_opto-std_opto), min(mean_flashopto-std_flashopto)]);

nexttile([1,1]);
patch('XData',[onset_i,onset_i,nbin,nbin],'YData',[minval,maxval,maxval,minval],...
    'FaceColor', [204/255,1,204/255], 'FaceAlpha', 0.5,'EdgeColor','none');
hold on;
stdshade_mean_std_linecolor(mean_flash,std_flash,0.25,[0.3,0.3,0.3],[0,0,0]);
ylim([minval,maxval]);
yticks([]);
xticks([]);
ax4=gca;
linkaxes([ax1,ax4],'x');

nexttile([1,1]);
patch('XData',[onset_i,onset_i,nbin,nbin],'YData',[minval,maxval,maxval,minval],...
    'FaceColor', [204/255,204/255,1], 'FaceAlpha', 0.5,'EdgeColor','none');
hold on;
stdshade_mean_std_linecolor(mean_opto,std_opto,0.25,[0.3,0.3,0.3],[0,0,0]);
ylim([minval,maxval]);
yticks([]);
xticks([]);
ax5=gca;
linkaxes([ax2,ax5],'x');

nexttile([1,1]);
patch('XData',[onset_i,onset_i,nbin,nbin],'YData',[minval,maxval,maxval,minval],...
    'FaceColor', [117, 192, 185]/255, 'FaceAlpha', 0.5,'EdgeColor','none');
hold on;
stdshade_mean_std_linecolor(mean_flashopto,std_flashopto,0.25,[0.3,0.3,0.3],[0,0,0]);
ylim([minval,maxval]);
yticks([]);
xticks([]);
ax6=gca;
linkaxes([ax3,ax6],'x');

n=size(flash_resp,1);
rec_id_u=unique(rec_id);
nrec=length(rec_id_u);
animal_i=zeros(1,nrec);
for i=1:nrec
    reci=rec_names(i);
    ii=find(contains([animal_ids(:).name], reci));
    animal_i(i)=animal_ids(ii).id;
end
animal_u=unique(animal_i);
nanimals=length(animal_u);
title(t,{titlestr; ['n=', num2str(n), '; nrec=', num2str(nrec), '; nanimals=', num2str(nanimals)]});

savefig(f,[figname,'.fig']);
saveas(f,[figname,'.svg']);
save([figname,'.mat'],'flash_resp','opto_resp','flashopto_resp','mean_flash','mean_opto','mean_flashopto',...
    'std_flash','std_opto','std_flashopto','bine','max_n','depth', 'nrec', "nanimals", 'rec_id','rec_names', 'animal_u');
end


function make_raster_all_responses(all_responses, recording_id, rec_names, animal_ids, bine, figname)
%% figure all neurons 

%normalize responses
max_n = max(all_responses(:,3:end),[],2);
all_responses_n=all_responses;
all_responses_n(:,3:end)=all_responses_n(:,3:end)./max_n;
%sort by depth
all_responses_n=sortrows(all_responses_n,2,"ascend");
maxd=all_responses_n(end,2);

nbin = length(bine)-1;
onset_i=find(bine>=0,1,"first");

f =figure; t=tiledlayout(1,6,"TileSpacing",'compact');
nexttile;
imagesc(all_responses_n(:,3:3+nbin-1), [0,1]); title('vis on');
nn=size(all_responses_n,1);
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);

depth_ticks=[0:100:maxd];
% for each depth increment find the index of the neuron
yticks_idx=arrayfun(@(i) find(all_responses_n(:,2)>=depth_ticks(i),1,"first"), 1:length(depth_ticks));
[yticks_idx,ia]=unique(yticks_idx);
t_ticks=1:10:nbin;
xticks_labels=bine(t_ticks);
yticks(yticks_idx);
yticklabels(depth_ticks(ia));
xticks(t_ticks);
xticklabels(xticks_labels);

nexttile;
imagesc(all_responses_n(:,3+nbin:3+2*nbin-1), [0,1]); title('vis off');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks(yticks_idx);
yticklabels({});
xticks(t_ticks);
xticklabels(xticks_labels);

nexttile;
imagesc(all_responses_n(:,3+2*nbin:3+3*nbin-1), [0,1]); title('vis+opto on');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks(yticks_idx);
yticklabels({});
xticks(t_ticks);
xticklabels(xticks_labels);

nexttile;
imagesc(all_responses_n(:,3+3*nbin:3+4*nbin-1), [0,1]); title('vis+opto off');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks(yticks_idx);
yticklabels({});
xticks(t_ticks);
xticklabels(xticks_labels);

nexttile;
imagesc(all_responses_n(:,3+4*nbin:3+5*nbin-1), [0,1]); title('opto on');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks(yticks_idx);
yticklabels({});
xticks(t_ticks);
xticklabels(xticks_labels);

nexttile;
imagesc(all_responses_n(:,3+5*nbin:end), [0,1]); title('opto off');
hold on;
plot([onset_i,onset_i],[0.5,nn+0.5],'--', 'Color',[64, 149, 206]/255);
yticks(yticks_idx);
yticklabels({});
xticks(t_ticks);
xticklabels(xticks_labels);

colormap(flipud(gray));
n=size(all_responses_n,1);
rec_id = unique(recording_id);
nrec=length(rec_id);
animal_i=zeros(1,nrec);
for i=1:nrec
    reci=rec_names(rec_id(i));
    ii=find(contains([animal_ids(:).name], reci));
    animal_i(i)=animal_ids(ii).id;
end
animal_u=unique(animal_i);
nanimals=length(animal_u);
title(t,{'all opto-flash reponsive units'; ['n=', num2str(n), '; nrec=', num2str(nrec), '; nanimals=', num2str(nanimals)]});

savefig(f,[figname,'.fig']);
saveas(f,[figname,'.svg']);
save([figname,'.mat'],'all_responses_n','bine', "nrec", "nanimals", 'recording_id','rec_names', 'animal_u');
end


function make_first_spike_timing_flash_opto_scatterplot(all_firstspike, recording_id, rec_names, animal_ids, minD, maxD, labelstr, path_to_figures)
% scatter plot of the spike count for flash off vs flash off+opto
% preselect neurons from the given depth range [minD,maxD)
% all_spikecount: clusterid, depth, flash_on, flash_off, flash_on+opto,
% flash_off+opto, opto on, opto off;


fs_flash = [];
fs_flashopto=[];
rec_id=[];

for i=1:size(all_firstspike,1)
    if ~(all_firstspike(i,2)>=minD && all_firstspike(i,2)<maxD)
        continue;
    end
    rec_id=[rec_id, recording_id(i)];
    fs_flash = [fs_flash, all_firstspike(i,4)];
    fs_flashopto=[fs_flashopto, all_firstspike(i,6)];
end

n= length(rec_id);
rec_id_u=unique(rec_id);
nrec=length(rec_id_u);
animal_i=zeros(1,nrec);
for i=1:nrec
    reci=rec_names(i);
    ii=find(contains([animal_ids(:).name], reci));
    animal_i(i)=animal_ids(ii).id;
end
animal_u=unique(animal_i);
nanimals=length(animal_u);
n_subtitle = ['n=', num2str(n), '; nrec=', num2str(nrec), '; nanimals=', num2str(nanimals)];


%test if the difference is significant
%wilcoxon signed rank test because pairs of first spike t come from the same
%neuron
alpha=0.01;
[p,h,stats]  = signrank(fs_flash,fs_flashopto,'tail','both','alpha',alpha);
signrank_testres.significance_level=alpha;
signrank_testres.tail='both'; 
signrank_testres.test_res_h=h;
signrank_testres.pvalue=p;
signrank_testres.zstatistics=stats.zval;
signrank_testres.stats=stats;

%make figure
f= figure;
mxt=max(max(fs_flash), max(fs_flashopto));
scatter(fs_flash,fs_flashopto,54,'filled', ...
    'MarkerFaceColor',[0.5, 0.5, 0.5],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on;
plot([0,mxt],[0,mxt],'--k','LineWidth',1);
axis equal;
ylim([0,mxt]);
xlim([0,mxt]);
xlabel('T first spike FLASH in [0, 0.2]s');
ylabel('T first spike FLASH+OPTO in [0, 0.2]s');
title({labelstr;...
['meandiff = ', num2str(mean(fs_flashopto-fs_flash,'omitnan')*1000), 'ms'];...
['mediandiff = ', num2str(median(fs_flashopto-fs_flash,'omitnan')*1000), 'ms'];...
['is diff: ', num2str(h), ' (zstat = ', num2str(round(stats.zval,3)),', p = ',num2str(p),')']; n_subtitle});


savefig(f,fullfile(path_to_figures,['T_first_spike_flash_opto_',labelstr,'.fig']));
saveas(f,fullfile(path_to_figures,['T_first_spike_flash_opto_',labelstr,'.svg']));
save(fullfile(path_to_figures,['T_first_spike_flash_opto_',labelstr,'_data.mat']),'fs_flash','fs_flashopto','signrank_testres', 'n', 'nrec', "nanimals", 'rec_id','rec_names', 'animal_u');

end


function make_spike_count_flash_opto_scatterplot(all_spikecount, recording_id, rec_names, animal_ids, minD, maxD, labelstr, path_to_figures)
% scatter plot of the spike count for flash off vs flash off+opto
% preselect neurons from the given depth range [minD,maxD)
% all_spikecount: clusterid, depth, flash_on, flash_off, flash_on+opto,
% flash_off+opto, opto on, opto off;


sc_flash = [];
sc_flashopto=[];
rec_id=[];

for i=1:size(all_spikecount,1)
    if ~(all_spikecount(i,2)>=minD && all_spikecount(i,2)<maxD)
        continue;
    end
    rec_id=[rec_id, recording_id(i)];
    sc_flash = [sc_flash, all_spikecount(i,4)];
    sc_flashopto=[sc_flashopto, all_spikecount(i,6)];
end

n= length(rec_id);
rec_id_u=unique(rec_id);
nrec=length(rec_id_u);
animal_i=zeros(1,nrec);
for i=1:nrec
    reci=rec_names(i);
    ii=find(contains([animal_ids(:).name], reci));
    animal_i(i)=animal_ids(ii).id;
end
animal_u=unique(animal_i);
nanimals=length(animal_u);
n_subtitle = ['n=', num2str(n), '; nrec=', num2str(nrec), '; nanimals=', num2str(nanimals)];


%test if the difference is significant
%wilcoxon signed rank test because pairs of observations are from one neuron
alpha=0.01;
[p,h,stats]  = signrank(sc_flash,sc_flashopto,'tail','both','alpha',alpha);
signrank_testres.significance_level=alpha;
signrank_testres.tail='both'; 
signrank_testres.test_res_h=h;
signrank_testres.pvalue=p;
signrank_testres.zstatistics=stats.zval;
signrank_testres.stats=stats;

f= figure;
mxt=max(max(sc_flash), max(sc_flashopto));
scatter(sc_flash,sc_flashopto,54,'filled', ...
    'MarkerFaceColor',[0.5, 0.5, 0.5],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on;
plot([0,mxt],[0,mxt],'--k','LineWidth',1);
axis equal;
ylim([0,mxt]);
xlim([0,mxt]);
xlabel('spike count FLASH [0, 0.2]s');
ylabel('spike count FLASH+OPTO [0, 0.2]s');
title({labelstr;...
['is diff: ', num2str(h), ' (zstat = ', num2str(round(stats.zval,3)),', p = ',num2str(p),')']; n_subtitle});

savefig(f,fullfile(path_to_figures,['spike_count_flash_opto_',labelstr,'.fig']));
saveas(f,fullfile(path_to_figures,['spike_count_flash_opto_',labelstr,'.svg']));
save(fullfile(path_to_figures,['spike_count_flash_opto_',labelstr,'_data.mat']),'sc_flash','sc_flashopto','signrank_testres', 'n', 'nrec', "nanimals", 'rec_id','rec_names', 'animal_u');
end


function make_spike_count_bl_blopto_scatterplot(all_spikecount, recording_id, rec_names, animal_ids, minD, maxD, labelstr, path_to_figures)
% scatter plot of the spike count for baseline vs flash baseline+opto
% preselect neurons from the given depth range [minD,maxD)
% all_spikecount: clusterid, depth, flash_on, flash_off, flash_on+opto,
% flash_off+opto, opto on, opto off;


sc_bl = [];
sc_blopto=[];
rec_id=[];

for i=1:size(all_spikecount,1)
    if ~(all_spikecount(i,2)>=minD && all_spikecount(i,2)<maxD)
        continue;
    end
    %% during 0.2s of opto
%     sc_bl = [sc_bl, all_spikecount(i,end-2)];
%     sc_blopto=[sc_blopto, all_spikecount(i,end-4)];
    %% during 0.1s of opto
    rec_id=[rec_id, recording_id(i)];
    sc_bl = [sc_bl, all_spikecount(i,end-1)];
    sc_blopto=[sc_blopto, all_spikecount(i,end)];
end

n= length(rec_id);
rec_id_u=unique(rec_id);
nrec=length(rec_id_u);
animal_i=zeros(1,nrec);
for i=1:nrec
    reci=rec_names(i);
    ii=find(contains([animal_ids(:).name], reci));
    animal_i(i)=animal_ids(ii).id;
end
animal_u=unique(animal_i);
nanimals=length(animal_u);
n_subtitle = ['n=', num2str(n), '; nrec=', num2str(nrec), '; nanimals=', num2str(nanimals)];


%test if the difference is significant
%wilcoxon signed rank test because pairs of observations are from one neuron
alpha=0.01;
[p,h,stats]  = signrank(sc_bl,sc_blopto,'tail','both','alpha',alpha);
signrank_testres.significance_level=alpha;
signrank_testres.tail='both'; 
signrank_testres.test_res_h=h;
signrank_testres.pvalue=p;
signrank_testres.zstatistics=stats.zval;
signrank_testres.stats=stats;

f= figure;
mxt=max(max(sc_bl), max(sc_blopto));
scatter(sc_bl,sc_blopto,54,'filled', ...
    'MarkerFaceColor',[0.5, 0.5, 0.5],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on;
plot([0,mxt],[0,mxt],'--k','LineWidth',1);
axis equal;
ylim([0,mxt]);
xlim([0,mxt]);
xlabel('spike count Baseline [0, 0.1]s');
ylabel('spike count Baseline+OPTO [0, 0.1]s');
title({labelstr;...
['is diff: ', num2str(h), ' (zstat = ', num2str(round(stats.zval,3)),', p = ',num2str(p),')']; n_subtitle});

savefig(f,fullfile(path_to_figures,['spike_count_bl_opto_',labelstr,'.fig']));
saveas(f,fullfile(path_to_figures,['spike_count_bl_opto_',labelstr,'.svg']));
save(fullfile(path_to_figures,['spike_count_bl_opto_',labelstr,'_data.mat']),'sc_bl','sc_blopto','signrank_testres', 'n', 'nrec', "nanimals", 'rec_id','rec_names', 'animal_u');
end






