%this script finds aborted exits as the trajectories near the border of the
%platform with subsequent reversal of the speed direction

% saveresfolder = 'D:\data\tnt_data\traj';
% load(fullfile(saveresfolder, 'allpts.mat'));

% save_res_folder = 'D:\data\tnt_data\traj\all';
% load(fullfile(save_res_folder, 'new_allpts.mat'));

% save_res_folder = 'D:\data\tnt_data\traj\all\vs_sham_only';
% load(fullfile('D:\data\tnt_data\traj\all', 'new_allpts.mat'));
clear all;
close all;
% save_res_folder = 'D:\data\tnt_data\traj\may2024_10min';
save_res_folder ='D:\data\resub\data\cliff\10min';
load(fullfile('D:\data\tnt_data\traj\all', 'new_allpts_updated.mat'));

crossings=allpts;
ta_s=2;
fps=60;
ta=round(ta_s*fps);

%max duration of recording to consider for the analysis
maxTmin=10;
maxT=maxTmin*60*fps;
kmax=Inf;
%if a crossing was detected skip ahead to avoid double-counting
fskip_after=round(1*fps);

%arena parameters
plst=cliffparams.area_size(1:2);
plen=cliffparams.area_size(3:4);
bdpos=cliffparams.platform_border; 
transfer_zone=cliffparams.platform_offset; 

%do we want plots per recording?
indiv_plots=true;
color_crossings=false; %set true if you want colored trajectories

aborted_tnt=[];
aborted_sham=[];
aborted_traj_tnt=[];
aborted_traj_sham=[];

crossing_speed={};

if color_crossings % set color indices to color the exit trajectories
    indcol=-ta:ta;
    col_traj=indcol/length(indcol);
end

skip_wt=true;
figfile=fullfile(save_res_folder,['aborted_traj_',num2str(maxTmin),'min.pdf']);
for i=1:length(crossings)
    if skip_wt && contains(crossings(i).animal_name,'wt','IgnoreCase',true)
        continue;
    end
    animal=crossings(i).animal_name;
    animal=animal(1:min(length(animal),50));
    animal_id=strrep(animal,'_',' ');

    rear_xy=crossings(i).xy(1:maxT,:); %xy of right ear of the mouse    
    
    %set the array for crossing traj
    aborted_traj=[];    

    cmap=turbo(201);
    lastT=min(maxT,size(rear_xy,1));
    kx=0;
    j=ta;
    jadd=false;
    if indiv_plots
        fig=figure;
    end
    % iterate thru all frames 
    while j<lastT-ta
        if j<ta || j>lastT-ta
            j=j+1;
            continue;
        end
        if kx>kmax
            break;
        end
        %position of the mouse on two consecutive frames
        x1=rear_xy(j,1);
        y1=rear_xy(j,2);
        x2=rear_xy(j+1,1);
        y2=rear_xy(j+1,2);

        %% check if platform border was crossed
        %crossing the top bd out: y1<bd & y2>=bd
        if y1<bdpos && y2>=bdpos
            x=(bdpos-y1)*(x2-x1)/(y2-y1)+x1; %where exactly crossing occured
            if x>bdpos, j=j+1; continue, end
            %mean speed before and after crossing
            mean_speed=mean(diff(rear_xy(j:j+ta,2),1));  % only y component
            mean_speed_b=mean(diff(rear_xy(j-ta:j,2),1));  % only y component
            if mean_speed_b<0, j=j+1; continue, end %entering from outside

            %how many frames before crossing inside the platform border
            pts_inside=sum(rear_xy(j-ta:j-1,2)<transfer_zone & ...
                rear_xy(j-ta:j-1,1)<transfer_zone);
            %aborted trajectory: if was inside the border and reversed speed
            if pts_inside>ta/2 && ...
                    mean_speed<0
                aborted_traj=cat(3,aborted_traj,rear_xy(j-ta:j+ta,:));
                if indiv_plots
                    if color_crossings, hold on; scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col_traj,'.');
                    else, hold on; plot(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),'k');
                    end
                end
                j=j+fskip_after; %skip frames after crossing to avoid double counting
                jadd=true;
            end            
        end
        %crossing the rigth bd out: x1<bd & x2>=bd
        if x1<bdpos && x2>=bdpos
            y=(bdpos-x1)*(y2-y1)/(x2-x1)+y1; %where exactly crossing occured
            if y>bdpos, j=j+1; continue, end
            %mean speed before and after crossing
            mean_speed_b=mean(diff(rear_xy(j-ta:j,1),1));  % only x component
            mean_speed=mean(diff(rear_xy(j:j+ta,1),1));  % only x component

            if mean_speed_b<0, j=j+1; continue, end %entering from outside

            %how many frames before crossing inside the platform border
            pts_inside=sum(rear_xy(j-ta:j-1,2)<transfer_zone & ...
                rear_xy(j-ta:j-1,1)<transfer_zone);
            %aborted trajectory: if was inside the border and reversed speed
            if pts_inside>ta/2 && ...
                    mean_speed<0

                aborted_traj=cat(3,aborted_traj,rear_xy(j-ta:j+ta,:));
                if indiv_plots
                    if color_crossings, hold on; scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col_traj,'.');
                    else, hold on; plot(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),'k');
                    end
                end
                %                 hold on; plot([plst(1),bdpos],[bdpos,bdpos],'--k');
                %                 hold on; plot([bdpos, bdpos],[plst(2),bdpos],'--k');
                %                 xlim([plst(1),plen(1)]); ylim([plst(2),plen(2)]);
                %                 title(j);
                j=j+fskip_after; %skip frames after crossing to avoid double counting
                jadd=true;
            end
        end
        
        if jadd % if crossing detected, framecounter was already advanced
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

    if contains(animal,'tnt','IgnoreCase',true)
        aborted_tnt=[aborted_tnt,size(aborted_traj,3)];
        aborted_traj_tnt=cat(3,aborted_traj_tnt, aborted_traj);
    else
        aborted_sham=[aborted_sham, size(aborted_traj,3)];
        aborted_traj_sham=cat(3,aborted_traj_sham, aborted_traj);
    end

end

%% significance tests on the the number of aborted exits
[aborted_test_sham_tnt.ranksum.pvalue, aborted_test_sham_tnt.ranksum.h0, ...
    aborted_test_sham_tnt.ranksum.stat] = ranksum(aborted_sham,aborted_tnt,'tail','right');
aborted_test_sham_tnt.ranksum.test_type = 'ranksum tail right';

%% plot aborted exits
fig=figure;
bar([1:2],[mean(aborted_sham),mean(aborted_tnt)],'FaceColor','none');
hold on;
xvals=[ones(1,numel(aborted_sham)),ones(1,numel(aborted_tnt))*2];
yvals=[aborted_sham,aborted_tnt];
s=swarmchart(xvals,yvals,'k','filled');
s.XJitterWidth=0.3;
xticks([1,2]);
xticklabels({'sham', 'tnt'});
title({'aborted crossings'; ...
    ['h0 =',num2str(aborted_test_sham_tnt.ranksum.h0), '; pvalue=',num2str(aborted_test_sham_tnt.ranksum.pvalue)]});
figfile=fullfile(save_res_folder,'sham_tnt_test_aborted_crossings');
exportgraphics(fig,[figfile,'.pdf'], 'BackgroundColor','none','ContentType','vector','Append', true);

matfile=fullfile(save_res_folder,'sham_tnt_test_aborted_crossings_data.mat');
save(matfile,'aborted_sham', 'aborted_tnt', 'aborted_traj_tnt','aborted_traj_sham','aborted_test_sham_tnt');

