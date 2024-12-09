function quantify_response_width_spikecount(pop1hist, pop2hist, depth1, depth2, pvals1, pvals2, pop1_ani, pop2_ani, pop1_recid, pop2_recid, histc, minD, maxD, pvalthresh, figsavename)
%plot spikecount and width of the responses of the two populations
% the width is estimated using the width of the fitted gaussian

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

%substract mean before the start of the stimulus and normalize max to 1
bin0=find(histc<0,1,'last');
mean_pop1_n=mean(pop1hist(:,1:bin0),2);
pop1hist_n=pop1hist-mean_pop1_n;
pop1hist_n=pop1hist_n./max(pop1hist_n,[],2);

mean_pop2_n=mean(pop2hist(:,1:bin0),2);
pop2hist_n=pop2hist-mean_pop2_n;
pop2hist_n=pop2hist_n./max(pop2hist_n,[],2);


%for each neuron with the max response after the sim onset fit the gaussian and
%collect the width, keep ids of ani and rec
[~,bin50]=min(abs(histc-0.05));
startt=histc(bin50);
starts=0.02; %0.04 ms response window
gaussEqn = 'exp(-(x-m)^2/(2*s^2))';
pop1_w=[];
pop1_selani_w=[];
pop1_selrec_w=[];
for ni=1:size(pop1hist_n,1)
    tri=pop1hist_n(ni,:);
    [~,tp]=max(tri);
    if tp<bin0, continue; end
%     tri(1:bin0)=0;
    startPoints = [startt,starts];
    [fgauss,gof,output] = fit(histc',tri',gaussEqn,'Start', startPoints,'Lower', [0,0], 'Upper', [0.4, 0.25]);
    pop1_w=[pop1_w,fgauss.s*2];
    pop1_selani_w=[pop1_selani_w, pop1_ani(ni)];
    pop1_selrec_w=[pop1_selrec_w, pop1_recid(ni)];
%     figure, plot(histc,tri); hold on; plot(fgauss);
end
pop2_w=[];
pop2_selani_w=[];
pop2_selrec_w=[];
for ni=1:size(pop2hist_n,1)
    tri=pop2hist_n(ni,:);
    [~,tp]=max(tri);
    if tp<bin0, continue; end
%     tri(1:bin0)=0;
    startPoints = [startt,starts];
    [fgauss,gof,output] = fit(histc',tri',gaussEqn,'Start', startPoints,'Lower', [0,0], 'Upper', [0.4, 0.25]);
    pop2_w=[pop2_w,fgauss.s*2];
    pop2_selani_w=[pop2_selani_w, pop2_ani(ni)];
    pop2_selrec_w=[pop2_selrec_w, pop2_recid(ni)];
%     figure, plot(histc,tri); hold on; plot(fgauss);
end

npop1_ani_w=numel(unique(pop1_selani_w));
npop2_ani_w=numel(unique(pop2_selani_w));
npop1_rec_w=numel(unique(pop1_selrec_w));
npop2_rec_w=numel(unique(pop2_selrec_w));
fig = figure;
% xvals=[ones(1,length(pop1_w)), ones(1,length(pop2_w))*2];
% yvals=[pop1_w, pop2_w];
% s=swarmchart(xvals, yvals,100,'filled');
s1=swarmchart(ones(1,length(pop1_w)), pop1_w,100,'filled', 'k'); hold on;
s2=swarmchart(ones(1,length(pop2_w))*2, pop2_w,100,'filled','r');
s1.XJitterWidth = 0.5; s1.MarkerFaceAlpha=0.2;
s2.XJitterWidth = 0.5; s2.MarkerFaceAlpha=0.2;
xticks([1,2]);
xticklabels({'Ctrl','TNT'});
hold on; plot(1,median(pop1_w),'+k','LineWidth',2);
plot(2,median(pop2_w),'+k','LineWidth',2);
ylabel('width of responses(s)');
[resp_diff_kstest.h,  resp_diff_kstest.p,  resp_diff_kstest.ksstat]=kstest2(pop1_w,pop2_w);
 resp_diff_kstest.name='kstest2';
titlestr={['kstest2: pvalue=', num2str( resp_diff_kstest.p),'; statistic: ', num2str( resp_diff_kstest.ksstat)];...
    ['minD=',num2str(minD),', maxD=',num2str(maxD), ', zetaThresh=',num2str(pvalthresh)];...
    ['ctrl: nn=',num2str(numel(pop1_w)), ', nrec=',num2str(npop1_rec_w), ', nani=',num2str(npop1_ani_w)];...
    ['tnt: nn=',num2str(numel(pop2_w)), ', nrec=',num2str(npop2_rec_w), ', nani=',num2str(npop2_ani_w)]};
title(titlestr);
exportgraphics(fig,figsavename, 'BackgroundColor','none','ContentType','vector', 'Append', true);

matsavename=[figsavename(1:end-4),'_width.mat'];
ctrl_resp_width=pop1_w;
tnt_resp_width=pop2_w;
save(matsavename, 'ctrl_resp_width','tnt_resp_width', 'resp_diff_kstest');


%alternatevely count number of spikes within 200ms window after stim onset for
%the neurons with the peak within this window
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
[ resp_diff_spikecount_kstest.h,  resp_diff_spikecount_kstest.p,  resp_diff_spikecount_kstest.ksstat]=kstest2(pop1_spc,pop2_spc);
 resp_diff_spikecount_kstest.name='kstest2';
titlestr={['kstest2: pvalue=', num2str( resp_diff_spikecount_kstest.p),'; statistic: ', num2str( resp_diff_spikecount_kstest.ksstat)];...
    ['minD=',num2str(minD),', maxD=',num2str(maxD), ', zetaThresh=',num2str(pvalthresh)];...
    ['ctrl: nn=',num2str(numel(pop1_spc)), ', nrec=',num2str(npop1_rec_spc), ', nani=',num2str(npop1_ani_spc)];...
    ['tnt: nn=',num2str(numel(pop2_spc)), ', nrec=',num2str(npop2_rec_spc), ', nani=',num2str(npop2_ani_spc)]};
title(titlestr);
exportgraphics(fig,figsavename, 'BackgroundColor','none','ContentType','vector', 'Append', true);

matsavename=[figsavename(1:end-4),'_spikecount.mat'];
ctrl_resp_spikecount=pop1_spc;
tnt_resp_spikecount=pop2_spc;
save(matsavename, 'ctrl_resp_spikecount','tnt_resp_spikecount', 'resp_diff_spikecount_kstest');
end