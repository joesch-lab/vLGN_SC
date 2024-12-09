%this script computes the time spent on the platform vs the time on the
%cliff; the strips near the walls as well as near the border of the
%platform are excluded from the analysis


clear all;
% close all;
save_res_folder = 'D:\data\tnt_data\traj\may2024_10min';
load(fullfile('D:\data\tnt_data\traj\all', 'new_allpts_updated.mat'));
%file with areana parameters: size of the areana and the platform;
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
bdpl_offset=[arenabd,arenabd];
bdpl_negoffset=[plbd,plbd];
bdrem_frac=cliffparams.remove_traj_bd_fraction; %0.1

sham_inds=[];
indiv_plots=false;
skip_wt=true;

for i=1:length(crossings)
    if skip_wt && contains(crossings(i).animal_name,'wt','IgnoreCase',true)
        continue;
    end  
    animal=crossings(i).animal_name;
    animal_id=strrep(animal,'_',' ');
    if contains(animal,'tnt','IgnoreCase',true)
        sham_inds=[sham_inds,0];
    else
        sham_inds=[sham_inds,1];
    end

    bdpl=[cliffparams.platform_border, cliffparams.platform_border];
    area_size=cliffparams.area_size(3)-cliffparams.area_size(1);
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

    %frames on the platform
    oni = sum(x>edvals(2) & x<bdpl_negoffset(2) & y>edvals(1) & y<bdpl_negoffset(1));
    %frames on the cliff
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

sham_inds=logical(sham_inds);
%% platform preference index
prefind = (onplatform_n'-offplatform_n')./(onplatform_n'+offplatform_n');
%% significance test
[p_rs, h0, ~] = ranksum(prefind(sham_inds),prefind(~sham_inds),"tail","right");

%% plot platform preference for tnt and sham
fig=figure;
bar([1:3],[mean(prefind(sham_inds)),NaN, mean(prefind(~sham_inds))],'FaceColor','none');
hold on;
xvals=[ones(1,sum(sham_inds)), 3*ones(1,sum(~sham_inds))];
s=swarmchart(xvals, [prefind(sham_inds)',prefind(~sham_inds)'],'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'sham','tnt'});
title({'platform preference index'; ...
    ['h0 =',num2str(h0), '; pvalue=',num2str(p_rs)]});
figfile=fullfile(save_res_folder,'sham_tnt_plpref');
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);

platform_preference.sham.onplatform=onplatform(sham_inds);
platform_preference.sham.offplatform=offplatform(sham_inds);
platform_preference.tnt.onplatform=onplatform(~sham_inds);
platform_preference.tnt.offplatform=offplatform(~sham_inds);

platform_preference.sham.onplatform_normed=onplatform_n(sham_inds);
platform_preference.sham.offplatform_normed=offplatform_n(sham_inds);
platform_preference.sham.prefind=prefind(sham_inds);
platform_preference.tnt.onplatform_normed=onplatform_n(~sham_inds);
platform_preference.tnt.offplatform_normed=offplatform_n(~sham_inds);
platform_preference.tnt.prefind=prefind(~sham_inds);

platform_preference.ranksum_test_prefind_diff.h0=h0;
platform_preference.ranksum_test_prefind_diff.pvalue=p_rs;
platform_preference.ranksum_test_prefind_diff.testtype ='ranksum tail right';

matfile=fullfile(save_res_folder,'sham_tnt_plpref.mat');
save(matfile,"platform_preference");


