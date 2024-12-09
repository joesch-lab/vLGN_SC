function mean_population_plots(pop1hist, pop2hist, depth1, depth2, pvals1, pvals2,  pop1_ani, pop2_ani, pop1_recid, pop2_recid, histc, minD, maxD, pvalthresh, centeronmax_yn, titlestr, figsavename)
% plot the mean responses of two populations

%if centered on max response is requested
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

fig =figure;
stdshade_mean_std_linecolor(mean_pop1,stderr_pop1,0.4,[0.5,0.5,0.5],[0,0,0],histc);
hold on;
stdshade_mean_std_linecolor(mean_pop2,stderr_pop2,0.4,[1,0,0],[1,0,0],histc);
if ~centeronmax_yn
    minval=min(min(mean_pop1-stderr_pop1),min(mean_pop2-stderr_pop2));
    minval=minval-0.05*minval;
    maxval=max(max(mean_pop1+stderr_pop1),max(mean_pop2+stderr_pop2));
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
exportgraphics(fig,figsavename, 'BackgroundColor','none','ContentType','vector', 'Append', true);
matsavename=[figsavename(1:end-4),'_meanresp.mat'];
save(matsavename, 'mean_pop1','stderr_pop1','mean_pop2','stderr_pop2','histc', 'pop1hist_n', 'pop2hist_n',... 
'pop1_ani', 'pop2_ani','pop1_recid', 'pop2_recid', 'minD', 'maxD', 'pvalthresh');
end