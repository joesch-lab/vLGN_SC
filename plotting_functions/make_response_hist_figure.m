function   make_response_hist_figure(spike_hist, histc, neuron_zetavals, depthSC, zetathresh, minD, maxD, titlestr, figsavename, sortby_t2maxFR, append_yn)
if ~exist('sortby_t2maxFR','var') || isempty(sortby_t2maxFR)
    sortby_t2maxFR=false;
end

if ~exist('zetathresh','var') || isempty(zetathresh)
    zetathresh=0.05;
end
if ~exist('append_yn','var') || isempty(append_yn)
    append_yn=true;
end

if isfield(neuron_zetavals,'z_pval')
    allpvals=[neuron_zetavals(:).z_pval];
else
    allpvals=neuron_zetavals(:);
end

%select neurons from specific depth
idx=find(depthSC>=minD & depthSC<maxD);
spike_hist=spike_hist(idx,:);
depthSC=depthSC(idx);
pvals_all=allpvals(idx);

%select responsive neurons
maxval=max(spike_hist,[],2);
idx=find(pvals_all<=zetathresh &maxval>0);
subtitle=['fraction of responsive neurons: ',num2str(numel(idx)/size(spike_hist,1))];
spike_hist=spike_hist(idx,:);
depthSC=depthSC(idx);
pvals_all=pvals_all(idx);


%make figure, where the neurons are sorted by depth, on the left there is
%depth bar, on the right there is a p value bar

%for visualization, the neuronal responses need to be normalized
spike_hist_normed=spike_hist;
maxval=max(spike_hist,[],2);
spike_hist_normed=spike_hist_normed./maxval;
if sortby_t2maxFR
    [~, idt] = max(spike_hist_normed,[],2);
    idt=histc(idt);
    sort_vals = [1:size(spike_hist,1);idt]';
else %sort by depth
    sort_vals=[1:size(spike_hist,1);depthSC']';
end
sort_vals=sortrows(sort_vals,2);

spike_hist_normed= spike_hist_normed(sort_vals(:,1),:);
depthSC=depthSC(sort_vals(:,1));
pvals_all=pvals_all(sort_vals(:,1));

population_average=mean(spike_hist_normed, "omitnan");

fig=figure; t=tiledlayout(4,6);

nexttile([3,4]);
imagesc(spike_hist_normed);
hold on;
nbin = length(histc);
bin0=interp1(histc,1:nbin,0);
plot([bin0,bin0],[0.5, size(spike_hist_normed,1)+0.5],'--b');
xticks(0.5:4:nbin+0.5);
xticklabels(histc(1:4:end));
if size(spike_hist_normed,1)<100,daylabel=4;
else, daylabel=20; end;
yticks(1:daylabel:size(spike_hist_normed,1));
yticklabels(sort_vals(1:daylabel:end,2));
colorbar('westoutside');
colormap(flipud(gray));

nexttile([3,2]); %pvalue bar
pvals_all_disp=min(pvals_all,0.07);
barh(pvals_all_disp);

yticks(1:daylabel:size(spike_hist_normed,1));
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
title(t, cat(1,depthtitle, titlestr, subtitle));


[filepath,~,~]=fileparts(figsavename);
matfile=fullfile(filepath,[strrep(titlestr{1},' ','_'),'.mat']);
save(matfile, 'spike_hist_normed', 'histc','pvals_all', 'depthSC', 'population_average', 'minD','maxD');
exportgraphics(fig,figsavename, 'BackgroundColor','none','ContentType','vector', 'Append', append_yn);
end