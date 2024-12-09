% load('D:\data\tnt_data\traj\sham_tnt_traj_data.mat');
% find the exits and compute the speed around them
saveresfolder = 'D:\data\tnt_data\traj';
ta_s=2;
fps=60;
ta=round(ta_s*fps);
bd=600;
ebd=200;
maxT=54000;
kmax=Inf;
fskip_after=round(1*fps);

cross_out_all_tnt=[];
cross_out_all_sham=[];
cross_out_c_tnt=[];
cross_out_c_sham=[];
cross_out_b_tnt=[];
cross_out_b_sham=[];

cross_out_all_tnt_normal=[];
cross_out_all_sham_normal=[];
cross_out_c_tnt_normal=[];
cross_out_c_sham_normal=[];
cross_out_b_tnt_normal=[];
cross_out_b_sham_normal=[];

aborted_tnt=[];
aborted_sham=[];
aborted_traj_tnt=[];
aborted_traj_sham=[];

crossing_speed={};

for i=1:8%length(crosssings)   
   animal=crosssings(i).animal_name;
   animal_id=strrep(animal,'_',' ');

   rear_xy=crosssings(i).xy(1:maxT,:);       

    if i==9
    rear_xy(1:404,:)=[];    
    end

    pos_filt=movmedian(rear_xy,[15,15],1);

    %interpolate nan values in ear positions
    nfr=size(rear_xy,1);
    inds=1:nfr;
    naninds=isnan(rear_xy(:,1)); 

    allvals=interp1(inds(~naninds),rear_xy(~naninds,1),inds,"linear","extrap");
    rear_xy(:,1)=allvals;
    naninds=isnan(rear_xy(:,2));    
    allvals=interp1(inds(~naninds),rear_xy(~naninds,2),inds,"linear","extrap");
    rear_xy(:,2)=allvals;

% %     if ismember(i,flipped_platform)
% %         maxx=max(rear_xy(:,1));
% %         rear_xy(:,1)=maxx-rear_xy(:,1);
% %     end


%     figure, scatter(rear_xy(:,1),rear_xy(:,2),'.k');

    plst=min(rear_xy);
    plen=max(rear_xy);
    plsize=plen-plst;
    bdpos=mean(plst+plsize/2);

    %set the array for crossing traj 
    cross_out_all=[];
    cross_out_c=[];
    cross_out_b=[];
    cross_out_all_normal=[];
    cross_out_c_normal=[];
    cross_out_b_normal=[];
    aborted_traj=[];
    indcol=-ta:ta;
    col_traj=indcol/length(indcol);

    exiti=1;    
    cmap=turbo(201);
    lastT=min(maxT,size(rear_xy,1));    
    kx=0;
    j=ta;
    jadd=false;
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
        
%         %crossing the top bd out: y1<bd & y2>=bd
%         if y1<bd && y2>=bd
%             x=(bd-y1)*(x2-x1)/(y2-y1)+x1;
%             if x>bd, continue, end
%             cross_out_all=cat(3,cross_out_all,rear_xy(j-ta:j+ta,:));
%             cross_speed_normal=diff(rear_xy(j-ta:j+ta,2),1);  % only y component
%             cross_out_all_normal=cat(2,cross_out_all_normal,cross_speed_normal);
%             if x<ebd
%                 cross_out_b=cat(3,cross_out_b,rear_xy(j-ta:j+ta,:));
%                 cross_out_b_normal=cat(2,cross_out_b_normal,cross_speed_normal);
%             else
%                 cross_out_c=cat(3,cross_out_c,rear_xy(j-ta:j+ta,:));
%                 cross_out_c_normal=cat(2,cross_out_c_normal,cross_speed_normal);
%                 kx=kx+1;
%             end
% %             figure, scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col,'.');
% 
%         end
% 
%         %crossing the rigth bd out: x1<bd & x2>=bd
%         if x1<bd && x2>=bd
%             y=(bd-x1)*(y2-y1)/(x2-x1)+y1;
%             if y>bd, continue, end
%             cross_out_all=cat(3,cross_out_all,rear_xy(j-ta:j+ta,:));
%             cross_speed_normal=diff(rear_xy(j-ta:j+ta,1),1);  % only x component
%             cross_out_all_normal=cat(2,cross_out_all_normal,cross_speed_normal);
%             if y<ebd
%                 cross_out_b=cat(3,cross_out_b,rear_xy(j-ta:j+ta,:));
%                 cross_out_b_normal=cat(2,cross_out_b_normal,cross_speed_normal);
%             else
%                 cross_out_c=cat(3,cross_out_c,rear_xy(j-ta:j+ta,:));
%                 cross_out_c_normal=cat(2,cross_out_c_normal,cross_speed_normal);
%                 kx=kx+1;
%             end
% %             figure, scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col,'.');
%         end

         %crossing the top bd out: y1<bd & y2>=bd
        if y1<bdpos && y2>=bdpos
            x=(bdpos-y1)*(x2-x1)/(y2-y1)+x1;
            if x>bdpos, j=j+1; continue, end            
            %mean speed after crossing
            mean_speed=mean(diff(rear_xy(j:j+ta,2),1));  % only y component
            mean_speed_b=mean(diff(rear_xy(j-ta:j,2),1));  % only y component            
            if mean_speed_b<0, j=j+1; continue, end %entering from outside
            if mean_speed<0
                aborted_traj=cat(3,aborted_traj,rear_xy(j-ta:j+ta,:));                
                figure, scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col_traj,'.');
                hold on; plot([plst(1),bdpos],[bdpos,bdpos],'--k'); 
                hold on; plot([bdpos, bdpos],[plst(2),bdpos],'--k'); 
                xlim([plst(1),plen(1)]); ylim([plst(2),plen(2)]);
                title(j);
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
            if mean_speed<0
                aborted_traj=cat(3,aborted_traj,rear_xy(j-ta:j+ta,:));                
                figure, scatter(rear_xy(j-ta:j+ta,1),rear_xy(j-ta:j+ta,2),[],col_traj,'.');
                hold on; plot([plst(1),bdpos],[bdpos,bdpos],'--k'); 
                hold on; plot([bdpos, bdpos],[plst(2),bdpos],'--k'); 
                xlim([plst(1),plen(1)]); ylim([plst(2),plen(2)]);
                title(j);
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

%     cross_speeds_all = squeeze(sum(diff(cross_out_all,1).^2,2).^0.5)';
%     cross_speeds_c = squeeze(sum(diff(cross_out_c,1).^2,2).^0.5)';
%     cross_speeds_b = squeeze(sum(diff(cross_out_b,1).^2,2).^0.5)';


%     crossing_speed(i).animal=animal;
%     crossing_speed(i).speed_all=cross_speeds_all;
%     crossing_speed(i).speed_center=cross_speeds_c;
%     crossing_speed(i).speed_bd=cross_speeds_b;
%     crossing_speed(i).speed_all_normal=cross_out_all_normal;
%     crossing_speed(i).speed_center_normal=cross_out_c_normal;
%     crossing_speed(i).speed_bd_normal=cross_out_b_normal;




%     fig=figure; t=tiledlayout(1,3,"TileSpacing","compact");
% cross_speeds_all = squeeze(sum(diff(cross_out_all,1).^2,2).^0.5);
% maxval=quantile(cross_speeds_all(:),0.95);
% nexttile;
% imagesc(cross_speeds_all',[0,maxval]);
% title('speed of all crossings');
% 
% cross_speeds_c = squeeze(sum(diff(cross_out_c,1).^2,2).^0.5);
% nexttile;
% imagesc(cross_speeds_c',[0,maxval]);
% title('speed of center crossings');
% 
% cross_speeds_b = squeeze(sum(diff(cross_out_b,1).^2,2).^0.5);
% nexttile;
% imagesc(cross_speeds_b',[0,maxval]);
% title('speed of border crossings');
% 
% title(t, animal);


if contains(animal,'tnt','IgnoreCase',true)
    cross_out_all_tnt=[cross_out_all_tnt; cross_speeds_all];
    cross_out_c_tnt=[cross_out_c_tnt;cross_speeds_c];
    cross_out_b_tnt=[cross_out_b_tnt;cross_speeds_b];
    cross_out_all_tnt_normal=[cross_out_all_tnt_normal; cross_out_all_normal'];
    cross_out_c_tnt_normal=[cross_out_c_tnt_normal;cross_out_c_normal'];
    cross_out_b_tnt_normal=[cross_out_b_tnt_normal;cross_out_b_normal'];
    aborted_tnt=[aborted_tnt,size(aborted_traj,3)];    
    aborted_traj_tnt=cat(3,aborted_traj_tnt, aborted_traj);
else
    cross_out_all_sham=[cross_out_all_sham;cross_speeds_all];
    cross_out_c_sham=[cross_out_c_sham;cross_speeds_c];
    cross_out_b_sham=[cross_out_b_sham;cross_speeds_b];
    cross_out_all_sham_normal=[cross_out_all_sham_normal;cross_out_all_normal'];
    cross_out_c_sham_normal=[cross_out_c_sham_normal;cross_out_c_normal'];
    cross_out_b_sham_normal=[cross_out_b_sham_normal;cross_out_b_normal'];
    aborted_sham=[aborted_sham, size(aborted_traj,3)];
    aborted_traj_sham=cat(3,aborted_traj_sham, aborted_traj);
end

end

[h_ac,p_ac,ks2stat_ac] = kstest2(aborted_sham,aborted_tnt);
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
    ['h0 sham=',num2str(h_ac), '; pvalue=',num2str(p_ac)]});
figfile=fullfile(save_res_folder,'sham_tnt_test_aborted_crossings');
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);


fig=figure; t=tiledlayout(2,3,"TileSpacing","compact");
maxval=quantile([abs(cross_out_all_tnt_normal(:)); abs(cross_out_all_sham_normal(:))],0.95);

%sham
nexttile;
imagesc(cross_out_all_sham_normal,[-maxval,maxval]);
title('sham - all');

nexttile;
imagesc(cross_out_c_sham_normal,[-maxval,maxval]);
title('center');

nexttile;
imagesc(cross_out_b_sham_normal,[-maxval,maxval]);
title('border');

%tnt
nexttile;
imagesc(cross_out_all_tnt_normal,[-maxval,maxval]);
title('tnt - all');

nexttile;
imagesc(cross_out_c_tnt_normal,[-maxval,maxval]);
title('center');

nexttile;
imagesc(cross_out_b_tnt_normal,[-maxval,maxval]);
title('border');

if kmax<Inf
    figfile=fullfile(saveresfolder,['exits_normal_vel_tnt_sham_first_',num2str(kmax)]);    
    matfile=fullfile(saveresfolder,['exits_normal_vel_tnt_sham_first_',num2str(kmax),'.mat']);    
else
    figfile=fullfile(saveresfolder,'exits_normal_vel_tnt_sham_all');
    matfile=fullfile(saveresfolder,'exits_normal_vel_tnt_sham_all.mat');    
end
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);

save(matfile,'cross_out_all_tnt_normal', 'cross_out_c_tnt_normal','cross_out_b_tnt_normal', ...
    'cross_out_all_sham_normal','cross_out_c_sham_normal', 'cross_out_b_sham_normal');



figure, plot(mean(cross_out_all_tnt_normal)); hold on; plot(mean(cross_out_all_sham_normal))


fig=figure; t=tiledlayout(2,3,"TileSpacing","compact");
maxval=quantile([cross_out_all_tnt(:); cross_out_all_sham(:)],0.95);

%sham
nexttile;
imagesc(cross_out_all_sham,[0,maxval]);
title('sham - all');

nexttile;
imagesc(cross_out_c_sham,[0,maxval]);
title('center');

nexttile;
imagesc(cross_out_b_sham,[0,maxval]);
title('border');

%tnt
nexttile;
imagesc(cross_out_all_tnt,[0,maxval]);
title('tnt - all');

nexttile;
imagesc(cross_out_c_tnt,[0,maxval]);
title('center');

nexttile;
imagesc(cross_out_b_tnt,[0,maxval]);
title('border');

if kmax<Inf
    figfile=fullfile(saveresfolder,['exits_speeds_tnt_sham_first_',num2str(kmax)]);    
    matfile=fullfile(saveresfolder,['exits_speeds_tnt_sham_first_',num2str(kmax),'.mat']);    
else
    figfile=fullfile(saveresfolder,'exits_speeds_tnt_sham_all');
    matfile=fullfile(saveresfolder,'exits_speeds_tnt_sham_all.mat');    
end
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);

save(matfile,'cross_out_all_tnt', 'cross_out_c_tnt','cross_out_b_tnt', ...
    'cross_out_all_sham','cross_out_c_sham', 'cross_out_b_sham');


figure, plot(mean(cross_out_c_sham)); hold on; plot(mean(cross_out_c_tnt))


