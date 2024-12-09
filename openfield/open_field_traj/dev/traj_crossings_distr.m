% repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\TNT_Sham_Max';
save_res_folder = 'D:\data\tnt_data\traj';
% 
% fl=dir(fullfile(repofolder ,'*.csv'));
% reminds=[];
% for i=1:length(fl)
%     if startsWith(fl(i).name,'._')
%         reminds=[reminds,i];
%     end
% end
% fl(reminds)=[];
% 
% %padded platform border
% bd=600;
% flipped_platform=5:9;
% platformbins=0:0.1:2;
% platformbins_1s=0.1:0.1:1;
platform_vals=[100:100:bd];
nb=numel(platform_vals);
% nbins1s=numel(platformbins_1s);
% 
% crosssings={};
% 
for i=1:length(fl)
   filename=fullfile(repofolder, fl(i).name);
   [~,animal,~]=fileparts(fl(i).name);
    animal_id=strrep(animal,'_',' ');

   opts = detectImportOptions(filename);
    opts.SelectedVariableNames = {'Var5','Var6','Var7','Var17','Var18','Var20','Var21','Var23','Var24','Var26','Var27'};%'pupil_el','pupil_az'};
    posxy = readmatrix(filename, opts);
    posxy(1:2,:)=[];
    bottom_left_top_right=[median(posxy(:,9),'omitnan'), median(posxy(:,4),'omitnan'),...
        median(posxy(:,5),'omitnan'),median(posxy(:,6),'omitnan')];
    rear_xy=posxy(:,1:2);       
    rear_xy_ll=posxy(:,3);       
    
    %interpolate nan values in ear positions
    nfr=size(rear_xy,1);
    inds=1:nfr;
    naninds=isnan(rear_xy(:,1) | rear_xy_ll<0.8);    
    allvals=interp1(inds(~naninds),rear_xy(~naninds,1),inds,"linear","extrap");
    rear_xy(:,1)=allvals;
    naninds=isnan(rear_xy(:,2));    
    allvals=interp1(inds(~naninds),rear_xy(~naninds,2),inds,"linear","extrap");
    rear_xy(:,2)=allvals;

    if ismember(i,flipped_platform)
        maxx=max(rear_xy(:,1));
        rear_xy(:,1)=maxx-rear_xy(:,1);
    end

%     figure, scatter(rear_xy(:,1),rear_xy(:,2),'.k');

    %set the bins for counting crossings        
%     cross_out_counts=zeros(1, numel(platformbins)-1);
%     cross_in_counts=zeros(1, numel(platformbins)-1);
    cross_out_counts=zeros(1, nb*2);
    cross_in_counts=zeros(1, nb*2);
    exiti=1;
    fig = figure;
    cmap=turbo(201);

    for j=1:nfr-1
        x1=rear_xy(j,1); 
        y1=rear_xy(j,2);
        x2=rear_xy(j+1,1); 
        y2=rear_xy(j+1,2);  
        
        %crossing the top bd out: y1<bd & y2>=bd
        if y1<bd && y2>=bd 
            x=(bd-y1)*(x2-x1)/(y2-y1)+x1;
            if x>bd, continue, end            
            [~, binind]=find(platform_vals-x>=0,1,"first");
            cross_out_counts(binind)=cross_out_counts(binind)+1;

            xplot=rear_xy(j-100:j+100,1);           
            yplot=rear_xy(j-100:j+100,2);          

            hold on;
            plot(xplot,yplot,'k');
%             scatter(xplot,yplot,15,cmap,'filled');


%             t=((bd-y1)*(x2-x1)/(y2-y1)+x1)/bd; %where the crossing happened
%             if t>1, continue; end            
%             binval=ceil(t*10)/10; %find the bin index
%             if binval<=0
%                 binind=1;
%             else
%                 binind=find(platformbins_1s==binval);
%             end
%             cross_out_counts(binind)=cross_out_counts(binind)+1;
        end

        %crossing the rigth bd out: x1<bd & x2>=bd
        if x1<bd && x2>=bd  
            y=(bd-x1)*(y2-y1)/(x2-x1)+y1;
            if y>bd, continue, end
            [~, binind]=find(platform_vals-y>=0,1,"first");
            binind=nb-binind+1;
            cross_out_counts(nb+binind)=cross_out_counts(nb+binind)+1;
%             t=1-((bd-x1)*(y2-y1)/(x2-x1)+y1)/bd; %where the crossing happened            
%             if t<0, continue; end            
%             binval=ceil(t*10)/10; %find the bin index
%             if binval>=1
%                 binind=nbins1s;
%             else
%                 binind=find(platformbins_1s==binval);
%             end
%             cross_out_counts(nbins1s+binind)=cross_out_counts(nbins1s+binind)+1;
        hold on;
        xplot=rear_xy(j-100:j+100,1);           
        yplot=rear_xy(j-100:j+100,2);   
        plot(xplot,yplot,'k');
%         scatter(xplot,yplot,15,cmap,'filled');
         

        end

        %crossing the top bd in: y1>=bd & y2<bd
         if y1>=bd && y2<bd
             x=(bd-y1)*(x2-x1)/(y2-y1)+x1;
            if x>bd, continue, end            
            [~, binind]=find(platform_vals-x>=0,1,"first");
            cross_in_counts(binind)=cross_in_counts(binind)+1;            
%             t=((bd-y1)*(x2-x1)/(y2-y1)+x1)/bd; %where the crossing happened
%             if t>1, continue; end            
%             binval=ceil(t*10)/10; %find the bin index
%             if binval<=0
%                 binind=1;
%             else
%                 binind=find(platformbins_1s==binval);
%             end
%             cross_in_counts(binind)=cross_in_counts(binind)+1;
        end

        %crossing the rigth bd in: x1>=bd & x2<bd
        if x1>=bd && x2<bd   
            y=(bd-x1)*(y2-y1)/(x2-x1)+y1;
            if y>bd, continue, end            
            [~, binind]=find(platform_vals-y>=0,1,"first");
            binind=nb-binind+1;
            cross_in_counts(nb+binind)=cross_in_counts(nb+binind)+1;
%             t=1-((bd-x1)*(y2-y1)/(x2-x1)+y1)/bd; %where the crossing happened
%             if t<0, continue; end            
%             binval=ceil(t*10)/10; %find the bin index
%             if binval>=1
%                 binind=nbins1s;
%             else
%                 binind=find(platformbins_1s==binval);
%             end
%             cross_in_counts(nbins1s+binind)=cross_in_counts(nbins1s+binind)+1;
        end
    end
% %     figure, plot(cross_out_counts); hold on; plot(cross_in_counts);
% %     hold on;
% %     maxval=max([cross_out_counts,cross_in_counts]);
% %     plot([nbins1s,nbins1s],[0,maxval],'--k');
% %     xticks([]);
% %     xlabel('top edge             rigth edge');    
% %     legend({'out','in'});    
% %     title(animal_id);
%     crosssings(i).in=cross_in_counts;
%     crosssings(i).out=cross_out_counts;
%     crosssings(i).animal_name=animal_id;
%     crosssings(i).binval=platformbins;
%     crosssings(i).xy=rear_xy;

% colorbar;
% caxis([-100, 100]);
xlim(bottom_left_top_right([2,4])); 
ylim(bottom_left_top_right([1,3]));
       
title(strrep(['exits for ',fl(i).name],'_',' '));
figfile=fullfile(save_res_folder,['exits_',crosssings(i).animal_name]);
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);
end





sham_in=[];
sham_out=[];
tnt_in=[];
tnt_out=[];
for i=1:length(crosssings)
    fig = figure(Position=[100,100,800,400]); t=tiledlayout(1,2);
    nexttile;
    scatter(crosssings(i).xy(:,1),crosssings(i).xy(:,2),'.k');
    maxval=max(crosssings(i).xy(:));
    xlim([1,maxval]); ylim([1,maxval]);

    nexttile;
    plot(crosssings(i).out); hold on; plot(crosssings(i).in);
    hold on;
    maxval=max([crosssings(i).out, crosssings(i).in]);
    plot([nb,nb],[0,maxval],'--k');
    xticks([]);
    xlabel('top edge             rigth edge');    
    legend({'out','in'});    
    title(t,crosssings(i).animal_name);

    figfile=fullfile(save_res_folder,[strrep(crosssings(i).animal_name,' ','_'),'_traj']);
    saveas(fig,[figfile,'.fig']);
    saveas(fig,[figfile,'.png']);
    saveas(fig,[figfile,'.svg']);

    if contains(crosssings(i).animal_name,'tnt',IgnoreCase=true)
        tnt_in=[tnt_in;crosssings(i).in];
        tnt_out=[tnt_out;crosssings(i).out];
    else
         sham_in=[sham_in;crosssings(i).in];
         sham_out=[sham_out;crosssings(i).out];
    end
end



sham_in_n=sham_in./sum(sham_in,2);
sham_out_n=sham_out./sum(sham_out,2);
tnt_in_n=tnt_in./sum(tnt_in,2);
tnt_out_n=tnt_out./sum(tnt_out,2);

mean_sham_in=mean(sham_in_n);
mean_sham_out=mean(sham_out_n);
mean_tnt_in=mean(tnt_in_n);
mean_tnt_out=mean(tnt_out_n);
stderr_sham_in=std(sham_in_n)/(size(sham_in_n,1)^0.5);
stderr_sham_out=std(sham_out_n)/(size(sham_out_n,1)^0.5);
stderr_tnt_in=std(tnt_in_n)/(size(tnt_in_n,1)^0.5);
stderr_tnt_out=std(tnt_out_n)/(size(tnt_out_n,1)^0.5);

sham_all=cat(1,sham_out,sham_in);
tnt_all=cat(1,tnt_in,tnt_out);
sham_all_n=cat(1,sham_in_n,sham_out_n);
tnt_all_n=cat(1,tnt_in_n,tnt_out_n);

mean_sham_all=mean(sham_all_n);
mean_tnt_all=mean(tnt_all_n);
stderr_sham_all=std(sham_all_n)/(size(sham_all_n,1)^0.5);
stderr_tnt_all=std(mean_tnt_all)/(size(mean_tnt_all,1)^0.5);

%% crossings out and in
fig = figure(Position=[100,100,800,400]); t=tiledlayout(1,2);
nexttile;
stdshade_mean_std_linecolor(mean_sham_out,stderr_sham_out,0.4,[0.5,0.5,0.5],[0,0,0]);
hold on;
stdshade_mean_std_linecolor(mean_tnt_out,stderr_tnt_out,0.4,[0.5,0,0],[1,0,0]);
hold on;
maxval=max([mean_sham_out+stderr_sham_out, mean_tnt_out+stderr_tnt_out]);
plot([nb,nb],[0,maxval],'--k');
xticks([]);
xlabel('top edge             rigth edge');
legend({'sham','tnt'});
title('crossing out');

nexttile;
stdshade_mean_std_linecolor(mean_sham_in,stderr_sham_in,0.4,[0.5,0.5,0.5],[0,0,0]);
hold on;
stdshade_mean_std_linecolor(mean_tnt_in,stderr_tnt_in,0.4,[0.5,0,0],[1,0,0]);
hold on;
maxval=max([mean_sham_in+stderr_sham_in, mean_tnt_in+stderr_tnt_in]);
plot([nb,nb],[0,maxval],'--k');
xticks([]);
xlabel('top edge             rigth edge');
legend({'sham','tnt'});
title('crossing in');

figfile=fullfile(save_res_folder,['sham_tnt_traj_nn']);
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);

%% all crossings together
fig = figure;
stdshade_mean_std_linecolor(mean_sham_all,stderr_sham_all,0.4,[0.5,0.5,0.5],[0,0,0]);
hold on;
stdshade_mean_std_linecolor(mean_tnt_all,stderr_tnt_all,0.4,[0.5,0,0],[1,0,0]);
hold on;
maxval=max([mean_sham_all+stderr_sham_all, mean_tnt_all+stderr_tnt_all]);
plot([nb,nb],[0,maxval],'--k');
xticks([]);
xlabel('top edge             rigth edge');
legend({'sham','tnt'});
title('all crossings');
figfile=fullfile(save_res_folder,['sham_tnt_traj_nn_allcrossings']);
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);


%% make the statistical test that the number of crossings 
% at the borders are different from the crossings in the center
binborder=[1,1,0,0,0,0,0,0,0,0,1,1]>0;
sham_border_center=cat(2,sum(sham_all(:,binborder),2)/sum(binborder),sum(sham_all(:,~binborder),2)/sum(~binborder));
tnt_border_center=cat(2,sum(tnt_all(:,binborder),2)/sum(binborder),sum(tnt_all(:,~binborder),2)/sum(~binborder));
[h_sham,p_sham,ks2stat_sham] = kstest2(sham_border_center(:,1),sham_border_center(:,2));
[h_tnt,p_tnt,ks2stat_tnt] = kstest2(tnt_border_center(:,1),tnt_border_center(:,2));
fig=figure;
bar([1:5],[mean(sham_border_center),NaN, mean(tnt_border_center)],'FaceColor','none');
hold on;
xvals=[ones(1,size(sham_border_center,1)), 2*ones(1,size(sham_border_center,1)),...
    4*ones(1,size(tnt_border_center,1)),5*ones(1,size(tnt_border_center,1))];
yvals=[sham_border_center(:,1)',sham_border_center(:,2)',tnt_border_center(:,1)',tnt_border_center(:,2)'];
s=swarmchart(xvals,yvals,'k','filled');
s.XJitterWidth=0.3;
xticks([1,2,4,5]);
xticklabels({'sham bd', 'sham cen','tnt bd', 'tnt cen'});
title({'all crossings'; ...
    ['h0 sham=',num2str(h_sham), '; pvalue=',num2str(p_sham)]; ......
    ['h0 tnt=',num2str(h_tnt), '; pvalue=',num2str(p_tnt)]});
figfile=fullfile(save_res_folder,'sham_tnt_test_allcrossings');
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);


%% same test for crossings out
sham_border_center_out=cat(2,sum(sham_out(:,binborder),2)/sum(binborder),sum(sham_out(:,~binborder),2)/sum(~binborder));
tnt_border_center_out=cat(2,sum(tnt_out(:,binborder),2)/sum(binborder),sum(tnt_out(:,~binborder),2)/sum(~binborder));
[h_sham_out,p_sham_out,ks2stat_sham_out] = kstest2(sham_border_center_out(:,1),sham_border_center_out(:,2));
[h_tnt_out,p_tnt_out,ks2stat_tnt_out] = kstest2(tnt_border_center_out(:,1),tnt_border_center_out(:,2));
fig=figure;
bar([1:5],[mean(sham_border_center_out),NaN, mean(tnt_border_center_out)],'FaceColor','none');
hold on;
xvals=[ones(1,size(sham_border_center_out,1)), 2*ones(1,size(sham_border_center_out,1)),...
    4*ones(1,size(tnt_border_center_out,1)),5*ones(1,size(tnt_border_center_out,1))];
yvals=[sham_border_center_out(:,1)',sham_border_center_out(:,2)',tnt_border_center_out(:,1)',tnt_border_center_out(:,2)'];
s=swarmchart(xvals,yvals,'k','filled');
s.XJitterWidth=0.3;
xticks([1,2,4,5]);
xticklabels({'sham bd', 'sham cen','tnt bd', 'tnt cen'});
title({'OUT crossings'; ...
    ['h0 sham=',num2str(h_sham_out), '; pvalue=',num2str(p_sham_out)]; ......
    ['h0 tnt=',num2str(h_tnt_out), '; pvalue=',num2str(p_tnt_out)]});
figfile=fullfile(save_res_folder,'sham_tnt_test_out');
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);

%% same test for crossings in
sham_border_center_in=cat(2,sum(sham_in(:,binborder),2)/sum(binborder),sum(sham_in(:,~binborder),2)/sum(~binborder));
tnt_border_center_in=cat(2,sum(tnt_in(:,binborder),2)/sum(binborder),sum(tnt_in(:,~binborder),2)/sum(~binborder));
[h_sham_in,p_sham_in,ks2stat_sham_in] = kstest2(sham_border_center_in(:,1),sham_border_center_in(:,2));
[h_tnt_in,p_tnt_in,ks2stat_tnt_in] = kstest2(tnt_border_center_in(:,1),tnt_border_center_in(:,2));
fig=figure;
bar([1:5],[mean(sham_border_center_in),NaN, mean(tnt_border_center_in)],'FaceColor','none');
hold on;
xvals=[ones(1,size(sham_border_center_in,1)), 2*ones(1,size(sham_border_center_in,1)),...
    4*ones(1,size(tnt_border_center_in,1)),5*ones(1,size(tnt_border_center_in,1))];
yvals=[sham_border_center_in(:,1)',sham_border_center_in(:,2)',tnt_border_center_in(:,1)',tnt_border_center_in(:,2)'];
s=swarmchart(xvals,yvals,'k','filled');
s.XJitterWidth=0.3;
xticks([1,2,4,5]);
xticklabels({'sham bd', 'sham cen','tnt bd', 'tnt cen'});
title({'IN crossings'; ...
    ['h0 sham=',num2str(h_sham_in), '; pvalue=',num2str(p_sham_in)]; ......
    ['h0 tnt=',num2str(h_tnt_in), '; pvalue=',num2str(p_tnt_in)]});
figfile=fullfile(save_res_folder,'sham_tnt_test_in');
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);

crossings_all.sham.out=sham_out;
crossings_all.sham.out_n=sham_out_n;
crossings_all.sham.out_mean=mean_sham_out;
crossings_all.sham.out_stderr=stderr_sham_out;
crossings_all.sham.in=sham_in;
crossings_all.sham.in_n=sham_in_n;
crossings_all.sham.in_mean=mean_sham_in;
crossings_all.sham.in_stderr=stderr_sham_in;
crossings_all.sham.all=sham_all;
crossings_all.sham.all_n=sham_all_n;
crossings_all.sham.all_mean=mean_sham_all;
crossings_all.sham.all_stderr=stderr_sham_all;

crossings_all.tnt.out=tnt_out;
crossings_all.tnt.out_n=tnt_out_n;
crossings_all.tnt.out_mean=mean_tnt_out;
crossings_all.tnt.out_stderr=stderr_tnt_out;
crossings_all.tnt.in=tnt_in;
crossings_all.tnt.in_n=tnt_in_n;
crossings_all.tnt.in_mean=mean_tnt_in;
crossings_all.tnt.in_stderr=stderr_tnt_in;
crossings_all.tnt.all=tnt_all;
crossings_all.tnt.all_n=tnt_all_n;
crossings_all.tnt.all_mean=mean_tnt_all;
crossings_all.tnt.all_stderr=stderr_tnt_all;

crossings_all.binvals=[0:100:1200];
crossings_all.planformedge=nb;

crossings_all.tests.all_crossings.sham_border_center=sham_border_center;
crossings_all.tests.all_crossings.tnt_border_center=tnt_border_center;
crossings_all.tests.all_crossings.sham_h0=h_sham;
crossings_all.tests.all_crossings.sham_pvalue=p_sham;
crossings_all.tests.all_crossings.testname='kstest2';
crossings_all.tests.all_crossings.tnt_h0=h_tnt;
crossings_all.tests.all_crossings.tnt_pvalue=p_tnt;

crossings_all.tests.out.sham_border_center=sham_border_center_out;
crossings_all.tests.out.tnt_border_center=tnt_border_center_out;
crossings_all.tests.out.sham_h0=h_sham_out;
crossings_all.tests.out.sham_pvalue=p_sham_out;
crossings_all.tests.out.testname='kstest2';
crossings_all.tests.out.tnt_h0=h_tnt_out;
crossings_all.tests.out.tnt_pvalue=p_tnt_out;

crossings_all.tests.in.sham_border_center=sham_border_center_in;
crossings_all.tests.in.tnt_border_center=tnt_border_center_in;
crossings_all.tests.in.sham_h0=h_sham_in;
crossings_all.tests.in.sham_pvalue=p_sham_in;
crossings_all.tests.in.testname='kstest2';
crossings_all.tests.in.tnt_h0=h_tnt_in;
crossings_all.tests.in.tnt_pvalue=p_tnt_in;


save(fullfile(save_res_folder,'sham_tnt_traj_data.mat'),'crosssings','crossings_all');

