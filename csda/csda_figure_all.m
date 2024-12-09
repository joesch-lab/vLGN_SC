function csda_figure_all(csda_info,csda_res_folder, vis_only)
% plot current sourse dencity raw data and dynamics of max gradient
% for 50ms before and 200ms after the onset
if ~exist('vis_only',"var") || isempty(vis_only)
    vis_only=0;
end

t_b_plot_s=0.05;
t_a_plot_s=0.2;
sr=20000;
t_b_plot=t_b_plot_s*sr;
t_a_plot=t_a_plot_s*sr;

rb=bluered();
ncha=length(csda_info.channel_labels);
cha_spacing = 25; %25 um between channels
onset_edge=csda_info.flash_onset_samples;
% t_dur=size(csda_info.CSD_vis_on,2);
dt_sec=0.05;
dt=dt_sec*sr;
inflection_label=csda_info.inflection_cha_above;
inflection_index=find(csda_info.channel_labels==inflection_label);
st_depth=275-(inflection_index-1)*25;
en_depth=st_depth+(ncha-1)*cha_spacing;
labels=st_depth:25:en_depth;

st=max(1,onset_edge-t_b_plot);
en=min(onset_edge+t_a_plot,size(csda_info.CSD_vis_on,2));
onset_edge_plot=t_b_plot+1;
t_dur=(en+-st+1);

nancha=isnan(csda_info.CSD_vis_on(:,1));
nnz_cha=sum(~nancha);
labels=labels(~nancha);
CSD_vis_plot=csda_info.CSD_vis_on(~nancha,st:en);
if ~vis_only, CSD_opto_plot=csda_info.CSD_opto_on(~nancha,st:en); end

fig = figure; t=tiledlayout(2,2,'TileSpacing','compact');
nexttile;
maxval=max(abs(CSD_vis_plot(:)));
imagesc(CSD_vis_plot,[-maxval,maxval]);
colormap(rb);
hold on;
plot([onset_edge_plot,onset_edge_plot],[0.5,nnz_cha+0.5],'--k');
colorbar;
yticks(2:2:nnz_cha);
yticklabels(labels(2:2:nnz_cha));
xticks(1:dt:t_dur);
xticklabels(round((0:dt_sec:(t_dur/sr))-onset_edge_plot/sr,2));
title('vis ON');

nexttile;
if ~vis_only
    maxval=max(abs(CSD_opto_plot(:)));
    imagesc(CSD_opto_plot,[-maxval,maxval]);
    colormap(rb);
    hold on;
    plot([onset_edge_plot,onset_edge_plot],[0.5,nnz_cha+0.5],'--k');
    colorbar;
    yticks(2:2:nnz_cha);
    yticklabels(labels(2:2:nnz_cha));
    xticks(1:dt:t_dur);
    xticklabels(round((0:dt_sec:(t_dur/sr))-onset_edge_plot/sr,2));
    title('opto ON ');
end


nexttile;
maxval=max(csda_info.var_t_vis);
plot(csda_info.var_t_vis/maxval,'Color',[0,0.8,0]);
if ~vis_only
    hold on;
    maxval=max(csda_info.var_t_opto);
    plot(csda_info.var_t_opto/maxval,'Color',[0,0,0.8]);
end
hold on;
plot([onset_edge_plot,onset_edge_plot],[0 1.2],'--k');
xticks(1:dt:t_dur);
xticklabels(round((0:dt_sec:(t_dur/sr))-onset_edge_plot/sr,2));
yticks([0,1]);
yticklabels([0,1]);
title('Vis and Opto variance dynamics');


nexttile;
probe_cha=zeros(32,2);
notnan=~isnan(csda_info.envelope_vis);
max_vis_fold=max(csda_info.envelope_vis(notnan));
probe_cha(notnan,1)=csda_info.envelope_vis(notnan)/max_vis_fold;
patch(-probe_cha(:,1), 1:32,[204/255,1,204/255]);
if ~vis_only
    max_opto_fold=max(csda_info.envelope_opto(notnan));
    probe_cha(notnan,2)=csda_info.envelope_opto(notnan)/max_opto_fold;
    hold on; patch(probe_cha(:,2),1:32,[204/255,204/255,1]);
end
set(gca,'YDir','reverse');
xlim([-2;2]);
xticks([-1;1]);
if ~vis_only
    xticklabels([max_vis_fold,max_opto_fold]);
else
    xticklabels([max_vis_fold,max_vis_fold]);
end
yticklabels(0*25:4*25:31*25);
yticks(0:4:31);
yticklabels(0*25:4*25:31*25);
ip_id=find(csda_info.channel_labels==csda_info.inflection_cha_above);
hold on; plot([-2;2],[ip_id,ip_id],'--','Color',[0.2 0.2 0.2])
title('vis & opto var rate');


[binfolder,~,~]=fileparts(csda_info.binfile);
[~,binbase,~]=fileparts(binfolder);
title(t,[binbase(1:15),'_',binbase(end-10:end)],'Interpreter','none');
if exist('csda_res_folder','var') && ~isempty(csda_res_folder)
    figname=[binbase,'_csda_all.pdf'];
    figname=fullfile(csda_res_folder,figname);
    %     saveas(fig,figname);
    exportgraphics (fig,figname, 'BackgroundColor','none','ContentType','vector');
end
end


