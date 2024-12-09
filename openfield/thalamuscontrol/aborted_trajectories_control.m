%this script finds aborted exits as the trajectories near the border of the
%platform with subsequent reversal of the speed

% saveresfolder = 'D:\data\tnt_data\traj';
% load(fullfile(saveresfolder, 'allpts.mat'));

% save_res_folder = 'D:\data\tnt_data\traj\all';
% load(fullfile(save_res_folder, 'new_allpts.mat'));

% save_res_folder = 'D:\data\tnt_data\traj\all\vs_sham_only';
% load(fullfile('D:\data\tnt_data\traj\all', 'new_allpts.mat'));

clear all;
close all;
save_res_folder = 'D:\data\tnt_data\traj\may2024_10min';
sham_matfile=fullfile(save_res_folder,'sham_tnt_test_aborted_crossings_data.mat');
load(sham_matfile);

save_res_folder_medThal = 'D:\data\resub\data\cliff\medial_thalamus';
figfile=fullfile(save_res_folder_medThal, 'aborted_traj_lin.pdf');
load(fullfile(save_res_folder_medThal, 'new_allpts_updated.mat'));
load('cliffparams.mat');
crossings=allpts;

ta_s=2;
fps=60;
ta=round(ta_s*fps);
bd=600;
ebd=200;
maxT=36000;%54000;
kmax=Inf;
fskip_after=round(1*fps);

aborted_medThal=[];
aborted_traj_medThal=[];

indiv_plots=true;
skip_wt=true;

for i=1:length(crossings)
    if skip_wt && contains(crossings(i).animal_name,'wt','IgnoreCase',true)
        continue;
    end
    animal=crossings(i).animal_name;
    animal=animal(1:min(length(animal),50));
    animal_id=strrep(animal,'_',' ');

    rear_xy=crossings(i).xy(1:maxT,:);


    plst=min(rear_xy);
    plen=max(rear_xy);
    plsize=plen-plst;
    bdpos=mean(plst+plsize/2);
    transfer_zone=bdpos+100;

    %set the array for crossing traj
    aborted_traj=[];
    indcol=-ta:ta;
    col_traj=indcol/length(indcol);

    exiti=1;
    cmap=turbo(201);
    lastT=min(maxT,size(rear_xy,1));
    kx=0;
    j=ta;
    jadd=false;
    if indiv_plots
        fig=figure;
    end
    while j<lastT-ta
        if j<ta || j>lastT-ta
            j=j+1;
            continue;
        end
        if kx>kmax
            break;
        end
        x1=rear_xy(j,1);
        y1=rear_xy(j,2);
        x2=rear_xy(j+1,1);
        y2=rear_xy(j+1,2);

        %crossing the top bd out: y1<bd & y2>=bd
        if y1<bdpos && y2>=bdpos
            x=(bdpos-y1)*(x2-x1)/(y2-y1)+x1;
            if x>bdpos, j=j+1; continue, end
            %mean speed after crossing
            mean_speed=mean(diff(rear_xy(j:j+ta,2),1));  % only y component
            mean_speed_b=mean(diff(rear_xy(j-ta:j,2),1));  % only y component
            if mean_speed_b<0, j=j+1; continue, end %entering from outside

            pts_inside=sum(rear_xy(j-ta:j-1,2)<transfer_zone & ...
                rear_xy(j-ta:j-1,1)<transfer_zone);
            if pts_inside>ta/2 && ...
                    mean_speed<0
                aborted_traj=cat(3,aborted_traj,rear_xy(j-ta:j+ta,:));
                if indiv_plots
                     hold on; plot(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2), 'k');
%                     hold on; scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col_traj,'.');
                end
                j=j+fskip_after;
                jadd=true;
            end
            %crossing the rigth bd out: x1<bd & x2>=bd
        end
        if x1<bdpos && x2>=bdpos
            y=(bdpos-x1)*(y2-y1)/(x2-x1)+y1;
            if y>bdpos, j=j+1; continue, end
            mean_speed_b=mean(diff(rear_xy(j-ta:j,1),1));  % only y component
            mean_speed=mean(diff(rear_xy(j:j+ta,1),1));  % only y component

            if mean_speed_b<0, j=j+1; continue, end %entering from outside

            %                if mean_speed<0
            pts_inside=sum(rear_xy(j-ta:j-1,2)<transfer_zone & ...
                rear_xy(j-ta:j-1,1)<transfer_zone);
            if pts_inside>ta/2 && ...
                    mean_speed<0

                aborted_traj=cat(3,aborted_traj,rear_xy(j-ta:j+ta,:));
                if indiv_plots
                    hold on; plot(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2), 'k');
%                     hold on; scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col_traj,'.');
                end
                %                 hold on; plot([plst(1),bdpos],[bdpos,bdpos],'--k');
                %                 hold on; plot([bdpos, bdpos],[plst(2),bdpos],'--k');
                %                 xlim([plst(1),plen(1)]); ylim([plst(2),plen(2)]);
                %                 title(j);
                j=j+fskip_after;
                jadd=true;
            end
        end
        if jadd
            jadd=false;
        else
            j=j+1;
        end
    end

    if indiv_plots
        hold on; plot([plst(1),bdpos],[bdpos,bdpos],'--k');
        hold on; plot([bdpos, bdpos],[plst(2),bdpos],'--k');
        xlim([plst(1),plen(1)]); ylim([plst(2),plen(2)]);
        title(animal_id);
        exportgraphics(fig,figfile, 'BackgroundColor','none','ContentType','vector','Append', true); 
    end
   
    aborted_medThal=[aborted_medThal,size(aborted_traj,3)];
    aborted_traj_medThal=cat(3,aborted_traj_medThal, aborted_traj);
end


[aborted_test_sham_medThalctrl.ranksum.pvalue, aborted_test_sham_medThalctrl.ranksum.h0, ...
    aborted_test_sham_medThalctrl.ranksum.stat] = ranksum(aborted_sham,aborted_medThal,'tail','right');
aborted_test_sham_medThalctrl.ranksum.test_type = 'ranksum tail right sham vs medThalctrl';

fig=figure;
bar([1:2],[mean(aborted_sham),mean(aborted_medThal)],'FaceColor','none');
hold on;
xvals=[ones(1,numel(aborted_sham)),ones(1,numel(aborted_medThal))*2];
yvals=[aborted_sham,aborted_medThal];
s=swarmchart(xvals,yvals,'k','filled');
s.XJitterWidth=0.3;
xticks([1,2]);
xticklabels({'sham', 'medThal ctrl'});
title({'aborted crossings'; ...
    ['h0 =',num2str(aborted_test_sham_medThalctrl.ranksum.h0), '; pvalue=',num2str(aborted_test_sham_medThalctrl.ranksum.pvalue)]});
exportgraphics(fig,figfile, 'BackgroundColor','none','ContentType','vector','Append', true); 

matfile=fullfile(save_res_folder_medThal,'sham_medThal_test_aborted_crossings_data.mat');
% save(matfile,'aborted_sham', 'aborted_medThal', 'aborted_traj_medThal','aborted_traj_sham','aborted_test_sham_medThalctrl');

