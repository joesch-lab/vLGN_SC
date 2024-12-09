%this script computes the time spent on the platform vs the time in the
%cliff; the strips near the walls as well as near the border of the
%platform are excluded from the analysis

% save_res_folder = 'D:\data\tnt_data\traj\all_animals_15min';
% load(fullfile(save_res_folder, 'new_allpts.mat'));

% save_res_folder = 'D:\data\tnt_data\traj\all';
% load(fullfile(save_res_folder, 'new_allpts.mat'));
clear all;
close all;
save_res_folder = 'D:\data\tnt_data\traj\may2024_10min';
sham_matfile=fullfile(save_res_folder,'sham_tnt_plpref.mat');
load(sham_matfile);

save_res_folder_ctrl = 'D:\data\resub\data\cliff\medial_thalamus';
load(fullfile(save_res_folder_ctrl, 'new_allpts_updated.mat'));
load('cliffparams.mat');
crossings=allpts;
fps=60; % bahav video fps

onplatform = [];
offplatform=[];
onplatform_n=[];
offplatform_n=[];
maxT=36000;%54000;

plbd=cliffparams.platform_border-cliffparams.platform_padding(1);%;500;
arenabd=cliffparams.platform_border+cliffparams.platform_padding(2);%500;
% plbd=550;
% arenabd=700;

bdpl_offset=[arenabd,arenabd];
bdpl_negoffset=[plbd,plbd];
sham_inds=[];
bdrem_frac=cliffparams.remove_traj_bd_fraction; %0.1

indiv_plots=false;
skip_wt=true;

for i=1:length(crossings)
    if skip_wt && contains(crossings(i).animal_name,'wt','IgnoreCase',true)
        continue;
    end  
    animal=crossings(i).animal_name;
    animal_id=strrep(animal,'_',' ');
   
    bdpl=[cliffparams.platform_border, cliffparams.platform_border];%bottom_left_top_right(1:2)+pl_size;
    area_size=cliffparams.area_size(3)-cliffparams.area_size(1);%;(mean(bottom_left_top_right(3:4))-mean(bottom_left_top_right(1:2)));
    edvals=cliffparams.area_size+bdrem_frac*area_size*[1,1,-1,-1];

    maxt=min(maxT,size(crossings(i).xy,1));

    x=crossings(i).xy(1:maxt,1);
    y=crossings(i).xy(1:maxt,2);

    indon=find(x>edvals(2) & x<bdpl_negoffset(2) & y>edvals(1) & y<bdpl_negoffset(1));
    indoff=find((x>bdpl_offset(2) & x<edvals (4) & y>edvals(1) & y<edvals(3)) | ...%rigth half
        (x>edvals(2) & x<=bdpl_offset(2) & y>bdpl_offset(1) & y<edvals(3)));
    if indiv_plots
        fig = figure; scatter(x,y,1,'.k'); hold on;
        scatter(x(indon),y(indon),1,'.b'); hold on;
        scatter(x(indoff),y(indoff),1,'.r'); hold on;
        title([crossings(i).animal_name]);
        figfile=fullfile(save_res_folder,[strrep(crossings(i).animal_name,' ','_'),'_prefind']);
        saveas(fig,[figfile,'.png']);
    end


    oni = sum(x>edvals(2) & x<bdpl_negoffset(2) & y>edvals(1) & y<bdpl_negoffset(1));

    offi=sum(x>bdpl_offset(2) & x<edvals (4) & y>edvals(1) & y<edvals(3)) + ...%rigth half
        sum(x>edvals(2) & x<=bdpl_offset(2) & y>bdpl_offset(1) & y<edvals(3));  %top left quater

    onplatform = [onplatform,oni];
    offplatform=[offplatform,offi];
    area1=(plbd-edvals(1))*(plbd-edvals(2));
    area2=(edvals(4)-arenabd)*(edvals(3)-edvals(1))+(arenabd-edvals(1))*(edvals(3)-arenabd);
    area_all=(edvals(4)-edvals(2))*(edvals(3)-edvals(1));
    onplatform_n = [onplatform_n,oni/(area1/area_all)];
    offplatform_n=[offplatform_n,offi/(area2/area_all)];
end

prefind = (onplatform_n'-offplatform_n')./(onplatform_n'+offplatform_n')

%test between sham animals and weak vLGN animals
sham_prefinds=platform_preference.sham.prefind;
[p_rs, h0, ~] = ranksum(sham_prefinds, prefind,"tail","right");

fig=figure;
bar([1:3],[mean(sham_prefinds),NaN, mean(prefind)],'FaceColor','none');
hold on;
xvals=[ones(1,numel(sham_prefinds)), 3*ones(1,numel(prefind))];
s=swarmchart(xvals, [sham_prefinds',prefind'],'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'sham','medThalamus ctrl'});
title({'platform preference index'; ...
    ['h0 =',num2str(h0), '; pvalue=',num2str(p_rs)]});
figfile=fullfile(save_res_folder_ctrl,'sham_medThalamus_plpref.pdf');
exportgraphics(fig,figfile, 'BackgroundColor','none','ContentType','vector','Append', true);    

platform_preference.medThalamus.onplatform=onplatform;
platform_preference.medThalamus.offplatform=offplatform;
platform_preference.medThalamus.onplatform_normed=onplatform_n;
platform_preference.medThalamus.offplatform_normed=offplatform_n;
platform_preference.medThalamus.prefind=prefind;

platform_preference.medThalamus.ranksum_test_prefind_diff.h0=h0;
platform_preference.medThalamus.ranksum_test_prefind_diff.pvalue=p_rs;
platform_preference.medThalamus.ranksum_test_prefind_diff.testtype ='ranksum tail right sham vs zi control';
matfile=fullfile(save_res_folder_ctrl,'sham_medThalamus_plpref.mat');
save(matfile,"platform_preference");


