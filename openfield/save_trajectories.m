%%  plot trajectories of the animal for the given amount of time
T2plot_min=10; %min
fps=60;
T2plot=T2plot_min*60*fps; %in sec

% % save_res_folder = 'D:\data\tnt_data\traj\may2024_temp';
% save_res_folder = 'D:\data\resub\data\cliff\medial_thalamus';
% load(fullfile(save_res_folder, 'new_allpts_updated.mat'));
% filout=fullfile(save_res_folder,['all_traj_', num2str(T2plot_min),'min.pdf']);

save_res_folder = 'D:\data\resub\data\cliff\10min';
load(fullfile(save_res_folder, 'new_allpts_updated.mat'));
filout=fullfile(save_res_folder,['all_traj_', num2str(T2plot_min),'min.pdf']);

for k=1:length(allpts)
    if contains(allpts(k).animal_name,'wt','IgnoreCase',true)
        continue;
    end
    rear_xy=  allpts(k).xy;
    animal_id=allpts(k).animal_name;

    maxTplot_i=min(T2plot, size(rear_xy,1));
    indscol=[1:maxTplot_i]/maxTplot_i;

    fig =figure; scatter(rear_xy(1:maxTplot_i,1),rear_xy(1:maxTplot_i,2),[],indscol,'.');
    title(animal_id(1:min(length(animal_id),50)));
    exportgraphics(fig,filout, 'BackgroundColor','none','ContentType','vector','Append', true);

    fig =figure; plot(rear_xy(1:maxTplot_i,1),rear_xy(1:maxTplot_i,2),'k');
    title(animal_id(1:min(length(animal_id),50)));
    exportgraphics(fig,filout, 'BackgroundColor','none','ContentType','vector','Append', true);
end