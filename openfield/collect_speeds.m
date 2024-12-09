%this script computes computes the speed of an animal, constructs speed 
% distributions and compares the speeds on the platform

%the strips near the walls are excluded from the analysis

clear all;
% close all;

save_res_folder = 'D:\data\tnt_data\traj\may2024_10min';
pix_scale=50/1000; %cm per pixel
load(fullfile('D:\data\tnt_data\traj\all', 'new_allpts_updated.mat'));
load('cliffparams.mat');

crossings=allpts;
fps=60; % behav video fps
tfit=0.25*fps; 

maxT=36000;

plbd=cliffparams.platform_border-cliffparams.platform_padding(1);%;500;
arenabd=cliffparams.platform_border+cliffparams.platform_padding(2);%500;
%fraction around border to be excluded
bdrem_frac=cliffparams.remove_traj_bd_fraction; %0.1 

sham_inds=[];
indiv_plots=false;
skip_wt=true;

%% iterate thru all animals
k=1;
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
    
    maxt=min(maxT,size(crossings(i).xy,1)); 
    tvals=1:maxt;
    outvals=zeros(size(tvals), 'logical');
    %find the outliers-times that animal travels too fast btw frames
    crossvals=crossings(i).xy(1:maxt,:);
    crossvals(:,1)=movmedian(crossvals(:,1),5);
    crossvals(:,2)=movmedian(crossvals(:,2),5);    
    crossings(i).xy(1:maxt,:)=crossvals;

    %% for each time point compute the speed fitted within 0.25sec
    speedsi=zeros(1,maxt-tfit);
    for j=1:maxt-tfit
        p0=crossings(i).xy(j,:);
        pj=crossings(i).xy(j+1:j+tfit,:);
        dj=(sum((pj-p0).^2,2)).^0.5;        
        %fitline
        coefficients = polyfit((1:tfit)/fps, dj, 1);
        % Get the slope.
        speedsi(j) = coefficients(1)*pix_scale;
    end

    speed_data(k).ani=crossings(i).animal_name;
    speed_data(k).speeds=speedsi;
    speed_data(k).xy=crossvals;
    k=k+1;
end
save(fullfile(save_res_folder,'speeds.mat'), 'speed_data', '-v7.3');

%% speed distributions pooled into tnt/sham groups 
load(fullfile(save_res_folder,'speeds.mat'));
tnt_samples=[];
ctrl_samples=[];
tnt_medians=[];
ctrl_medians=[];
area_size=cliffparams.area_size(3)-cliffparams.area_size(1);
edvals=cliffparams.area_size+bdrem_frac*area_size*[1,1,-1,-1];
pl_offset=cliffparams.platform_offset;

%% find extreme values to remove outliers
Mspeed=-Inf;
minspeed=Inf;
for i=1:length(speed_data)
    Mspeed=max(Mspeed, quantile(speed_data(i).speeds, 0.99));
    minspeed=min(minspeed, quantile(speed_data(i).speeds, 0.01));
end
%% collect all speed values per group and median values per animal
for i=1:length(speed_data)
    spi=speed_data(i).speeds;
    spi=spi(spi>minspeed); spi=spi(spi<Mspeed); %remove outliers

    if contains(speed_data(i).ani, 'tnt', 'IgnoreCase',true)
        tnt_medians=[tnt_medians,median(spi)];
        tnt_samples=[tnt_samples, spi];
    else
        ctrl_medians=[ctrl_medians, median(spi)];
        ctrl_samples=[ctrl_samples, spi];
    end
end
nbins=100;
bine=linspace(minspeed,Mspeed,nbins);
binc=bine(1:end-1)+diff(bine(1:2))/2;
hc_tnt=histcounts(tnt_samples,bine,'Normalization','probability');
hc_ctrl=histcounts(ctrl_samples,bine,'Normalization','probability');
fig=figure;
bar(binc,hc_ctrl, 'FaceColor','k',FaceAlpha=0.3); hold on;
bar(binc,hc_tnt, 'FaceColor','r',FaceAlpha=0.3);
[h0, p0,stati0]= kstest2(tnt_samples,ctrl_samples);
title({['speed distributions'];['kstest2: pvalue=', num2str(p0),'; kstats=', num2str(stati0)]});
xlabel('cm/s');
exportgraphics(fig,speeds_figname, 'BackgroundColor','none','ContentType','vector','Append', true);

fig=figure;
bar(1,mean(ctrl_medians), 'EdgeColor','k', FaceAlpha=0); hold on;
bar(3,mean(tnt_medians), 'EdgeColor','r', FaceAlpha=0); 
scatter(1+rand(size(ctrl_medians))*0.4-0.2, ctrl_medians, 'k', 'filled');
scatter(3+rand(size(tnt_medians))*0.4-0.2, tnt_medians, 'r', 'filled');
xticks([1,3]);
xticklabels({'ctrl', 'tnt'});
[p1, h1, stati1]= ranksum(ctrl_medians,tnt_medians);
title({['median speeds per animal'];['ranksum: pvalue=', num2str(p1),'; stats=', num2str(stati1.ranksum)]});
exportgraphics(fig,speeds_figname, 'BackgroundColor','none','ContentType','vector','Append', true);

%% save median values and test data
speeddata_grouped.all.medians.ctrl=ctrl_medians;
speeddata_grouped.all.medians.tnt=tnt_medians;
speeddata_grouped.all.vals.ctrl=ctrl_samples;
speeddata_grouped.all.vals.tnt=tnt_samples;
%tests on the difference between full distributions of speeds
speeddata_grouped.tests.all_distr_diff.testname='kstest2';
speeddata_grouped.tests.all_distr_diff.p= p0;
speeddata_grouped.tests.all_distr_diff.statval= stati0;
%tests on the difference between median values per animal
speeddata_grouped.tests.all_mediandiff.testname='ranksum';
speeddata_grouped.tests.all_mediandiff.p= p1;
speeddata_grouped.tests.all_mediandiff.statval= stati1;
save(fullfile(save_res_folder,'speeds_grouped.mat'), 'speeddata_grouped', '-v7.3');



