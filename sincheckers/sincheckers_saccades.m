%% finds saccades during sin checkers stimulus
% plots PSTH, mean poupulation response and 
% neuron spike count after saccades

%folder where to save data
resfolder='D:\data\resub\data\sin_checkers';

%folder where recordings with csd analysis with saved depth
csd_folder = 'D:\data\tnt_data\csda\'; 

%run sin_checkers to get the file eye_vals_per_rep_all.mat
load(fullfile(resfolder, 'eye_vals_per_rep_all.mat'));
figsavename=fullfile(resfolder, ['sacc_sincheckers.pdf']);

binw=0.01; %bin width
tbefore=0.5; %window before saccade
tafter=1; % time after saccade

% vectors to hold all neurons
sacc_hist=[];
aniid=[];
recid=[];
zeta=[];
tntid=[];
depth=[];
zetathresh=1;
plot_fig_perrecoding=false;

%% loop thru recordings, detect saccades, compute PSTH
for fi=1:length(sindata)
    clear data;

    currdate=sindata(fi).date(1:11);
    datfile=sindata(fi).datafile;

    load(datfile);
    disp(datfile);

    %% if recording doesn't have CSDA depth per neuron,
    % find another stimulus from the same
    %recording session and transfer CSD depth values
    if ~isfield(data.spike_data, 'SC_depth')
        %get experiment path
        if isfield(data.exp_info, 'folder')
            ofi=data.exp_info.folder;
        else
            ofi=data.exp_info.fullpath;
        end

        %find csd-file
        [~,recname,~]=fileparts(ofi);
        recwords=strsplit(recname,{'_','-'});
        numr=cellfun(@numel, recwords);
        recwords=recwords(numr>1);

        %find best matching csd file
        maxmatch=0;
        matid=-1;
        csda_list=dir(fullfile(csd_folder,'*.mat'));
        for ci=1:length(csda_list)
            maxmatchi=0;
            for wi=1:length(recwords)
                if contains(csda_list(ci).name,recwords(wi))
                    maxmatchi=maxmatchi+1;
                end
            end
            if maxmatchi>=maxmatch
                maxmatch=maxmatchi;
                matid=ci;
            end
        end

        csdfile=fullfile(csda_list(matid).folder, csda_list(matid).name);
        %assign depth values
        load(csdfile);
        %set depth of the neuron using csd
        data.csda=csda_data;
        data  = set_csda_depth_per_cluster(data);
        %resave data
        save(datfile, 'data');
    end

    %get spikes from all units
    spikes=data.spike_data;
    if data.trigger_data.cam_trigger(2)>100
        camtrig_s=data.trigger_data.cam_trigger/data.trigger_data.sampling_rate;
    else
        camtrig_s=data.trigger_data.cam_trigger;
    end

    minfrnum = min(numel(data.trigger_data.cam_trigger), numel(data.dlc.pupil_az));

    %find saccades
    eye_az= data.dlc.pupil_az - median(data.dlc.pupil_az,'omitnan');
    eye_el= data.dlc.pupil_el - median(data.dlc.pupil_el,'omitnan');
    camtrig_s=camtrig_s(1:minfrnum);
    eye_az= eye_az(1:minfrnum);
    eye_el= eye_el(1:minfrnum);

    [saccade,direction,speed,selBehavTimes, params] = get_saccades([],  camtrig_s ,eye_az, eye_el,[],[]);
    nsacc=length(saccade.PeakTime);
    disp(['Number of saccades: ',num2str(nsacc)]);
    sacc_framenum = zeros(nsacc,1);
    for si=1:nsacc
        sacc_framenum(si)=numel(find(camtrig_s<=saccade.PeakTime(si)));
    end
    saccade.Framenum=sacc_framenum;

    %collect spikes around saccade times
    if isfield(data.trigger_data,'neuronal_sampling_rate')
        sr_spiking=data.trigger_data.neuronal_sampling_rate;
    else
        sr_spiking=20000;
    end
    [saccade_spike_hist, histc] = get_spikes_on_times(saccade.PeakTime,spikes, sr_spiking, tbefore, tafter, binw);

    %make zeta-test for the responsiveness to saccades
    neuron_zetavals = stim_resp_zetatest(saccade.PeakTime,spikes, sr_spiking, 0.4);

    %depth of neurons from CSDA
    depthSC = [spikes(:).SC_depth]';

    %vector with the recording, animal id, istnt
    nn=size(saccade_spike_hist,1);
    recvec=ones(nn,1)*fi;
    aniid_fi=ones(nn,1)*sindata(fi).ani_id;
    istntvec=ones(nn,1)*sindata(fi).istnt;

    sacc_hist=[sacc_hist;saccade_spike_hist];
    aniid=[aniid;aniid_fi];
    recid=[recid;recvec];
    tntid=[tntid;istntvec];
    zeta=[zeta;[neuron_zetavals(:).z_pval]'];
    depth=[depth;depthSC];

    if plot_fig_perrecoding
        %make title for the figure
        titlestr=['animal id: ',num2str(sindata(fi).ani_id), '; ',currdate];
        if sindata(fi).istnt, titlestr=[titlestr, '; TNT'];
        else, if sindata(fi).is_wt, titlestr=[titlestr, '; WT'];
        else, titlestr=[titlestr, '; sham']; end; end
           
        %make figure of neuronal responses to saccades
        make_saccade_response_figure(saccade_spike_hist, histc, neuron_zetavals, depthSC, recvec, aniid_fi, zetathresh, -Inf, Inf, titlestr, figsavename, true, true, false);    
    end
end
save(fullfile(resfolder, 'spikehist_saccades.mat'),'sacc_hist', 'aniid', 'recid','tntid', 'zeta', 'depth', 'histc');


%% plots
load(fullfile(resfolder, 'spikehist_saccades.mat'));
% resfolder='D:\data\resub\data\sin_checkers\saccadic_suppression';
%% PSTH around saccades, tnt and sham, sSC and dSC
%%select tnt
tnt_ani=aniid(tntid>0);
tnt_sacc_hist=sacc_hist(tntid>0,:);
tnt_recid=recid(tntid>0);
tnt_zeta=zeta(tntid>0);
tnt_depth=depth(tntid>0);
% %make title for the figure
zthresh=0.01;
titlestr0=['TNT, nani: ',num2str(numel(unique(tnt_ani))), ', nrec: ',num2str(numel(unique(tnt_recid))),...
    ', nn: ',num2str(sum(tnt_zeta<=zthresh))];
titlestr2 = ['saccadic suppression, zeta pvalue = ',num2str(zthresh)];
titlestr={titlestr0; titlestr2};%;titlestr2;titlestr3};
%make figure of neuronal responses to saccades
figsavename=fullfile(resfolder, ['sacc_sincheckers_tnt_sSC.pdf']);
make_saccade_response_figure(tnt_sacc_hist, histc, tnt_zeta, tnt_depth, tnt_recid, tnt_ani, zthresh, -Inf, 500, 'TeLC', figsavename, true, false);
figsavename=fullfile(resfolder, ['sacc_sincheckers_tnt_dSC.pdf']);
make_saccade_response_figure(tnt_sacc_hist, histc, tnt_zeta, tnt_depth, tnt_recid, tnt_ani, zthresh, 500, Inf, 'TeLC', figsavename, true, false);
%
%%select ctrl
ctrl_ani=aniid(tntid<1);
ctrl_sacc_hist=sacc_hist(tntid<1,:);
ctrl_recid=recid(tntid<1);
ctrl_zeta=zeta(tntid<1);
ctrl_depth=depth(tntid<1);
% %make title for the figure
zthresh=0.01;
titlestr0=['ctrl, nani: ',num2str(numel(unique(ctrl_ani))), ', nrec: ',num2str(numel(unique(ctrl_recid))),...
    ', nn: ',num2str(sum(ctrl_zeta<=zthresh))];
titlestr2 = ['saccadic suppression, zeta pvalue = ',num2str(zthresh)];
titlestr={titlestr0; titlestr2};%;titlestr2;titlestr3};
%make figure of neuronal responses to saccades
figsavename=fullfile(resfolder, ['sacc_sincheckers_ctrl_sSC.pdf']);
make_saccade_response_figure(ctrl_sacc_hist, histc, ctrl_zeta, ctrl_depth, ctrl_recid, ctrl_ani, zthresh, -Inf, 500, 'CTRL', figsavename, true,false);
figsavename=fullfile(resfolder, ['sacc_sincheckers_ctrl_dSC.pdf']);
make_saccade_response_figure(ctrl_sacc_hist, histc, ctrl_zeta, ctrl_depth, ctrl_recid, ctrl_ani, zthresh, 500, Inf, 'CTRL', figsavename, true, false);

%% make a plot of population means
titlestr='Ctrl/TNT mean responses to saccades: aligned at the peak';
figsavename=fullfile(resfolder, ['mean_sacc_resp_peakaligned.pdf']);
mean_population_plots(ctrl_sacc_hist, tnt_sacc_hist, ctrl_depth, tnt_depth, ctrl_zeta, tnt_zeta, ctrl_ani, tnt_ani, ctrl_recid, tnt_recid, histc, -Inf, 500, zthresh, true, titlestr, figsavename);
%% make a violin plot of spike count for each neuron
titlestr='Ctrl/TNT mean responses to saccades: aligned to saccades';
figsavename=fullfile(resfolder, ['mean_sacc_resp_saccaligned.pdf']);
mean_population_plots(ctrl_sacc_hist, tnt_sacc_hist, ctrl_depth, tnt_depth, ctrl_zeta, tnt_zeta, ctrl_ani, tnt_ani, ctrl_recid, tnt_recid, histc, -Inf, 500, zthresh, false, titlestr, figsavename);
% figsavename=fullfile(resfolder, ['sacc_resp_width.pdf']);
figsavename=fullfile(resfolder, ['sacc_resp_spikecount.pdf']);
quantify_response_width(ctrl_sacc_hist, tnt_sacc_hist, ctrl_depth, tnt_depth, ctrl_zeta, tnt_zeta, ctrl_ani, tnt_ani, ctrl_recid, tnt_recid, histc, -Inf, 500, zthresh, figsavename);

%make plot of mean sham and tnt responses
function quantify_response_width(pop1hist, pop2hist, depth1, depth2, pvals1, pvals2, pop1_ani, pop2_ani, pop1_recid, pop2_recid, histc, minD, maxD, pvalthresh, figsavename)

if ~exist('pvalthresh','var') || isempty(pvalthresh)
    pvalthresh=0.05;
end

%select responsive neurons
if isfield(pvals1,'z_pval')
    pvals1=[pvals1(:).z_pval];
end
if isfield(pvals2,'z_pval')
    pvals2=[pvals2(:).z_pval];
end

%select neurons by depth and zeta value
idk=depth1>=minD & depth1<maxD & pvals1 <= pvalthresh;
pop1hist=pop1hist(idk,:); pop1_ani=pop1_ani(idk); pop1_recid=pop1_recid(idk);
idk=depth2>=minD & depth2<maxD & pvals2 <= pvalthresh;
pop2hist=pop2hist(idk,:); pop2_ani=pop2_ani(idk); pop2_recid=pop2_recid(idk);

%% count number of spikes within 200ms window after saccade
bin0=find(histc<0,1,'last');
bin200=find(histc<=0.2,1,'last');
binw=diff(histc(1:2));
pop1_spc=[];
pop1_selani_spc=[];
pop1_selrec_spc=[];
for ni=1:size(pop1hist,1)
    tri=pop1hist(ni,:);
    spc_i=sum(tri(bin0:bin200))*binw; %convert from fr/s to spike count
    pop1_spc=[pop1_spc,spc_i];
    pop1_selani_spc=[pop1_selani_spc, pop1_ani(ni)];
    pop1_selrec_spc=[pop1_selrec_spc, pop1_recid(ni)];
end
pop2_spc=[];
pop2_selani_spc=[];
pop2_selrec_spc=[];
for ni=1:size(pop2hist,1)
    tri=pop2hist(ni,:);
    spc_i=sum(tri(bin0:bin200))*binw;%convert from fr/s to spike count
    pop2_spc=[pop2_spc,spc_i];
    pop2_selani_spc=[pop2_selani_spc, pop2_ani(ni)];
    pop2_selrec_spc=[pop2_selrec_spc, pop2_recid(ni)];
end
npop1_ani_spc=numel(unique(pop1_selani_spc));
npop2_ani_spc=numel(unique(pop2_selani_spc));
npop1_rec_spc=numel(unique(pop1_selrec_spc));
npop2_rec_spc=numel(unique(pop2_selrec_spc));

fig = figure;
s1=swarmchart(ones(1,length(pop1_spc)), pop1_spc,100,'filled', 'k'); hold on;
s2=swarmchart(ones(1,length(pop2_spc))*2, pop2_spc,100,'filled','r');
s1.XJitterWidth = 0.5; s1.MarkerFaceAlpha=0.2;
s2.XJitterWidth = 0.5; s2.MarkerFaceAlpha=0.2;
xticks([1,2]);
xticklabels({'Ctrl','TNT'});
hold on; plot(1,median(pop1_spc),'+k','LineWidth',2);
plot(2,median(pop2_spc),'+k','LineWidth',2);
ylabel('spike count within 200ms');
[sacc_resp_diff_spikecount_kstest.h, sacc_resp_diff_spikecount_kstest.p, sacc_resp_diff_spikecount_kstest.ksstat]=kstest2(pop1_spc,pop2_spc);
sacc_resp_diff_spikecount_kstest.name='kstest2';
titlestr={['kstest2: pvalue=', num2str(sacc_resp_diff_spikecount_kstest.p),'; statistic: ', num2str(sacc_resp_diff_spikecount_kstest.ksstat)];...
    ['minD=',num2str(minD),', maxD=',num2str(maxD), ', zetaThresh=',num2str(pvalthresh)];...
    ['ctrl: nn=',num2str(numel(pop1_spc)), ', nrec=',num2str(npop1_rec_spc), ', nani=',num2str(npop1_ani_spc)];...
    ['tnt: nn=',num2str(numel(pop2_spc)), ', nrec=',num2str(npop2_rec_spc), ', nani=',num2str(npop2_ani_spc)]};
title(titlestr);
figsavename_spc=[figsavename(1:end-4),'_spikecount.pdf'];
exportgraphics(fig,figsavename_spc, 'BackgroundColor','none','ContentType','vector');

matsavename=[figsavename(1:end-4),'_spikecount.mat'];
ctrl_resp_spikecount=pop1_spc;
tnt_resp_spikecount=pop2_spc;
save(matsavename, 'ctrl_resp_spikecount','tnt_resp_spikecount', 'sacc_resp_diff_spikecount_kstest');
end

%make plot of mean sham and tnt responses
function mean_population_plots(pop1hist, pop2hist, depth1, depth2, pvals1, pvals2,  pop1_ani, pop2_ani, pop1_recid, pop2_recid, histc, minD, maxD, pvalthresh, centeronmax_yn, titlestr, figsavename)

if ~exist('centeronmax_yn','var') || isempty(centeronmax_yn)
    centeronmax_yn=false;
end
if ~exist('pvalthresh','var') || isempty(pvalthresh)
    pvalthresh=0.05;
end

%select responsive neurons
if isfield(pvals1,'z_pval')
    pvals1=[pvals1(:).z_pval];
end
if isfield(pvals2,'z_pval')
    pvals2=[pvals2(:).z_pval];
end

%select neurons by depth
idk=depth1>=minD & depth1<maxD & pvals1 <= pvalthresh;
pop1hist=pop1hist(idk,:); pop1_ani=pop1_ani(idk); pop1_recid=pop1_recid(idk);
idk=depth2>=minD & depth2<maxD & pvals2 <= pvalthresh;
pop2hist=pop2hist(idk,:); pop2_ani=pop2_ani(idk); pop2_recid=pop2_recid(idk);

bin0=find(histc<0,1,'last');
mean_pop1_n=mean(pop1hist(:,1:bin0),2);
pop1hist_n=pop1hist-mean_pop1_n;
pop1hist_n=pop1hist_n./max(pop1hist_n,[],2);

mean_pop2_n=mean(pop2hist(:,1:bin0),2);
pop2hist_n=pop2hist-mean_pop2_n;
pop2hist_n=pop2hist_n./max(pop2hist_n,[],2);

%if centered on the peak, shift histograms
if centeronmax_yn
    [~,peakt]=max(pop1hist_n,[],2);
    dt_peak=bin0-peakt;
    for ri=1:size(pop1hist_n,1)
        pop1hist_n(ri,:)=circshift(pop1hist_n(ri,:),dt_peak(ri));
    end

    [~,peakt]=max(pop2hist_n,[],2);
    dt_peak=bin0-peakt;
    for ri=1:size(pop2hist_n,1)
        pop2hist_n(ri,:)=circshift(pop2hist_n(ri,:),dt_peak(ri));
    end
end

mean_pop1=mean(pop1hist_n,'omitnan');
mean_pop2=mean(pop2hist_n,'omitnan');

stderr_pop1=std(pop1hist_n, 'omitnan')/(size(pop1hist_n,1)^0.5);
stderr_pop2=std(pop2hist_n, 'omitnan')/(size(pop2hist_n,1)^0.5);

bin1=-0.1;
bin1id=find(histc>=bin1,1,"first");
bin2=0.3;
bin2id=find(histc<=bin2,1,"last");


fig =figure;
stdshade_mean_std_linecolor(mean_pop1(bin1id:bin2id),stderr_pop1(bin1id:bin2id),0.4,[0.5,0.5,0.5],[0,0,0],histc(bin1id:bin2id));
hold on;
stdshade_mean_std_linecolor(mean_pop2(bin1id:bin2id),stderr_pop2(bin1id:bin2id),0.4,[1,0,0],[1,0,0],histc(bin1id:bin2id));
if ~centeronmax_yn
    minval=min(min(mean_pop1(bin1id:bin2id)-stderr_pop1(bin1id:bin2id)),min(mean_pop2(bin1id:bin2id)-stderr_pop2(bin1id:bin2id)));
    minval=minval-0.05*minval;
    maxval=max(max(mean_pop1(bin1id:bin2id)+stderr_pop1(bin1id:bin2id)),max(mean_pop2(bin1id:bin2id)+stderr_pop2(bin1id:bin2id)));
    maxval=maxval+0.05*maxval;
    plot([0,0], [minval,maxval],'--k');
end

npop1_ani=numel(unique(pop1_ani));
npop2_ani=numel(unique(pop2_ani));
npop1_rec=numel(unique(pop1_recid));
npop2_rec=numel(unique(pop2_recid));

titlestr={titlestr;...
    ['minD=',num2str(minD),', maxD=',num2str(maxD), ', zetaThresh=',num2str(pvalthresh)];...
    ['ctrl: nn=',num2str(size(pop1hist_n,1)), ', nrec=',num2str(npop1_rec), ', nani=',num2str(npop1_ani)];...
    ['tnt: nn=',num2str(size(pop2hist_n,1)), ', nrec=',num2str(npop2_rec), ', nani=',num2str(npop2_ani)]};
title(titlestr);
exportgraphics(fig,figsavename, 'BackgroundColor','none','ContentType','vector');
matsavename=[figsavename(1:end-3),'mat'];
save(matsavename, 'mean_pop1','stderr_pop1','mean_pop2','stderr_pop2','histc', 'pop1hist_n', 'pop2hist_n',... 
'pop1_ani', 'pop2_ani','pop1_recid', 'pop2_recid', 'minD', 'maxD', 'pvalthresh');
end


function   make_saccade_response_figure(saccade_spike_hist, histc, neuron_zetavals, depthSC, reclist, anilist, zetathresh, minD, maxD, titlename, figsavename, sortby_t2maxFR, append_yn, savedata)
if ~exist('sortby_t2maxFR','var') || isempty(sortby_t2maxFR)
    sortby_t2maxFR=false;
end

if ~exist('zetathresh','var') || isempty(zetathresh)
    zetathresh=0.05;
end
if ~exist('append_yn','var') || isempty(append_yn)
    append_yn=true;
end
if ~exist('savedata','var') || isempty(savedata)
    savedata=true;
end


if isfield(neuron_zetavals,'z_pval')
    allpvals=[neuron_zetavals(:).z_pval];
else
    allpvals=neuron_zetavals(:);
end
if size(allpvals,2)>size(allpvals,1)
    allpvals=allpvals';
end

%select neurons from specific depth
idx=find(depthSC>=minD & depthSC<maxD);
saccade_spike_hist=saccade_spike_hist(idx,:);
depthSC=depthSC(idx);
pvals_all=allpvals(idx);
reclist = reclist(idx);
anilist = anilist(idx);

%select responsive neurons
maxval=max(saccade_spike_hist,[],2);
idx=find(pvals_all<=zetathresh & maxval>0);
subtitle=['fraction of responsive neurons: ',num2str(numel(idx)/size(saccade_spike_hist,1))];
saccade_spike_hist=saccade_spike_hist(idx,:);
depthSC=depthSC(idx);
pvals_all=pvals_all(idx);
reclist = reclist(idx);
anilist = anilist(idx);


%make figure, where the neurons are sorted by depth, on the left there is
%depth bar, on the right there is a p value bar

%for visualization, the neuronal responses need to be normalized
saccade_spike_hist_normed=saccade_spike_hist;
maxval=max(saccade_spike_hist,[],2);
saccade_spike_hist_normed=saccade_spike_hist_normed./maxval;
if sortby_t2maxFR
    [~, idt] = max(saccade_spike_hist_normed,[],2);
    idt=histc(idt);
    sort_vals = [1:size(saccade_spike_hist,1);idt]';
else %sort by depth
    sort_vals=[1:size(saccade_spike_hist,1);depthSC']';
end
sort_vals=sortrows(sort_vals,2);

saccade_spike_hist_normed= saccade_spike_hist_normed(sort_vals(:,1),:);
bin1=-0.1;
bin1id=find(histc>=bin1,1,"first");
bin2=0.3;
bin2id=find(histc<=bin2,1,"last");
depthSC=depthSC(sort_vals(:,1));
pvals_all=pvals_all(sort_vals(:,1));

population_average=mean(saccade_spike_hist_normed, "omitnan");

fig=figure; t=tiledlayout(4,6);
population_average=population_average(bin1id:bin2id);
saccade_spike_hist_normed=saccade_spike_hist_normed(:,bin1id:bin2id);
saccade_spike_hist=saccade_spike_hist(:,bin1id:bin2id);
histc=histc(bin1id:bin2id);

nexttile([3,4]);
imagesc(saccade_spike_hist_normed);
hold on;
nbin = length(histc);
bin0=interp1(histc,1:nbin,0);
plot([bin0,bin0],[0.5, size(saccade_spike_hist_normed,1)+0.5],'--w');
xticks(0.5:4:nbin+0.5);
xticklabels(histc(1:4:end));
if size(saccade_spike_hist_normed,1)<100,daylabel=4;
else, daylabel=20; end;
yticks(1:daylabel:size(saccade_spike_hist_normed,1));
yticklabels(sort_vals(1:daylabel:end,2));
colormap(flipud(gray));
colorbar('westoutside');

nexttile([3,2]); %pvalue bar
pvals_all_disp=min(pvals_all,0.07);
barh(pvals_all_disp);

yticks(1:daylabel:size(saccade_spike_hist_normed,1));
yticklabels(sort_vals(1:daylabel:end,2));
hold on;
plot([0.05,0.05],[1,length(pvals_all)]);
set(gca,'YDir','reverse');
title('p-values');

nexttile([1,4]);%population plot
plot(histc,population_average); hold on;
minvp=min(population_average);
maxvp=max(population_average);
plot([0,0],[minvp-0.1*minvp, maxvp+0.1*maxvp],'--k');
xticks(histc(1:4:end));
ylim([[minvp-0.1*minvp, maxvp+0.1*maxvp]]);

depthtitle=['Depth: [',num2str(minD), ',',num2str(maxD),']'];
titlestr0=[titlename, ': nani= ',num2str(numel(unique(anilist))), ', nrec=',num2str(numel(unique(reclist))),...
     ', nn=',num2str(size(saccade_spike_hist_normed,1))];
titlestr2 = ['saccadic suppression, zeta pvalue = ',num2str(zetathresh)];
titlestr={titlestr0; titlestr2};

title(t, cat(1,depthtitle, titlestr, subtitle));
exportgraphics(fig,figsavename, 'BackgroundColor','none','ContentType','vector', 'Append', append_yn);
if savedata
    matsavename=[figsavename(1:end-4),'.mat'];
    save(matsavename, 'saccade_spike_hist_normed', 'saccade_spike_hist','histc',...
    'population_average','depthSC', 'pvals_all', 'reclist', 'anilist', '-v7.3');
end
end