
onplatform = [];
offplatform=[];
onplatform_n=[];
offplatform_n=[];
maxt=54000;

for i=1:length(crosssings)
    filename=fullfile(repofolder, fl(i).name);
    [~,animal,~]=fileparts(fl(i).name);
    animal_id=strrep(animal,'_',' ');

    opts = detectImportOptions(filename);
    opts.SelectedVariableNames = {'Var5','Var6','Var7','Var17','Var18','Var20','Var21','Var23','Var24','Var26','Var27'};%'pupil_el','pupil_az'};
    posxy = readmatrix(filename, opts);
    posxy(1:2,:)=[];
    bottom_left_top_right=[median(posxy(:,9),'omitnan'), median(posxy(:,4),'omitnan'),...
        median(posxy(:,5),'omitnan'),median(posxy(:,6),'omitnan')];

    pl_size=(mean(bottom_left_top_right(3:4))-mean(bottom_left_top_right(1:2)))/2;
    bdpl=bottom_left_top_right(1:2)+pl_size;    

    area_size=(mean(bottom_left_top_right(3:4))-mean(bottom_left_top_right(1:2)));
%     bdpl_offset=bottom_left_top_right(1:2)+pl_size+0.1*area_size;
    bdpl_offset=[600,600];
    bdpl_negoffset=[500,500];
%     bdpl_negoffset=bottom_left_top_right(1:2)+pl_size-0.1*area_size;
    edvals=[bottom_left_top_right]+0.1*area_size*[1,1,-1,-1];

    x=crosssings(i).xy(1:maxt,1);
    y=crosssings(i).xy(1:maxt,2);

    indon=find(x>edvals(2) & x<bdpl_negoffset(2) & y>edvals(1) & y<bdpl_negoffset(1));
    indoff=find((x>bdpl_offset(2) & x<edvals (4) & y>edvals(1) & y<edvals(3)) | ...%rigth half
        (x>edvals(2) & x<=bdpl_offset(2) & y>bdpl_offset(1) & y<edvals(3)));
    fig = figure; scatter(x,y,'.k'); hold on;
    scatter(x(indon),y(indon),'.b'); hold on;
    scatter(x(indoff),y(indoff),'.r'); hold on;
    title([crosssings(i).animal_name]);
    figfile=fullfile(save_res_folder,[strrep(crosssings(i).animal_name,' ','_'),'_prefind']);
    saveas(fig,[figfile,'.png']);


%     oni=sum(x<bdpl(2) & y<bdpl(1));
%     oni = sum(x>edvals(2) & x<bdpl(2) & y>edvals(1) & y<bdpl(1));
  oni = sum(x>edvals(2) & x<bdpl_negoffset(2) & y>edvals(1) & y<bdpl_negoffset(1));

    offi=sum(x>bdpl_offset(2) & x<edvals (4) & y>edvals(1) & y<edvals(3)) + ...%rigth half
        sum(x>edvals(2) & x<=bdpl_offset(2) & y>bdpl_offset(1) & y<edvals(3));  %top left quater 
%     nt=length(x)

    onplatform = [onplatform,oni];
    offplatform=[offplatform,offi];   
    area1=(500-edvals(1))*(500-edvals(2));
    area1=(600-edvals(1))*(600-edvals(2));
    area2=(edvals(4)-600)*(edvals(3)-edvals(1))+(600-edvals(1))*(edvals(3)-600);
    area_all=(edvals(4)-edvals(2))*(edvals(3)-edvals(1));
    onplatform_n = [onplatform_n,oni/(area1/area_all)];
    offplatform_n=[offplatform_n,offi/(area2/area_all)];   

%     figure, scatter(x,y,'.k');
end

% onplatform_n=onplatform/0.09;
% offplatform_n=offplatform/0.36;
% 
% % onplatform_n=onplatform/0.25;
% % % onplatform_n=onplatform/0.09;
% % offplatform_n=offplatform/0.75;
% on_off_n=[onplatform_n;offplatform_n]'
% on_off=[onplatform;offplatform]'

 prefind = (onplatform_n'-offplatform_n')./(onplatform_n'+offplatform_n')
[h0, pval,~] = kstest2(prefind(1:4),prefind(5:end));

fig=figure;
bar([1:3],[mean(prefind(1:4)),NaN, mean(prefind(5:end))],'FaceColor','none');
hold on;
xvals=[ones(1,4), 3*ones(1,5)];    
s=swarmchart(xvals,prefind,'k','filled');
s.XJitterWidth=0.3;
xticks([1,3]);
xticklabels({'sham','tnt'});
title({'platform preference index'; ...
    ['h0 =',num2str(h0), '; pvalue=',num2str(pval)]});
figfile=fullfile(save_res_folder,'sham_tnt_plpref');
saveas(fig,[figfile,'.fig']);
saveas(fig,[figfile,'.png']);
saveas(fig,[figfile,'.svg']);




% for i=4%1:8%length(crosssings)
%     h2=round(size(crosssings(i).xy,1)/2);
%     maxval=max(crosssings(i).xy(:));
%     x1=crosssings(i).xy(1:h2,1);
%     y1=crosssings(i).xy(1:h2,2);
%     x2=crosssings(i).xy(h2+1:end,1);
%     y2=crosssings(i).xy(h2+1:end,2);
%     
%     fig1 = figure;   
%     plot(x1,y1,'k');    
%     xlim([1,maxval]); ylim([1,maxval]);   
%     title([crosssings(i).animal_name, ' part1']);
%     figfile=fullfile(save_res_folder,[strrep(crosssings(i).animal_name,' ','_'),'_linetraj_1']);
% %     saveas(fig,[figfile,'.fig']);
%     saveas(fig1,[figfile,'.png']);
%     saveas(fig1,[figfile,'.svg']);
%     saveas(fig1,[figfile,'.pdf']);
% 
% 
%     fig2 = figure;   
%     plot(x2,y2,'k');    
%     xlim([1,maxval]); ylim([1,maxval]);   
%     title([crosssings(i).animal_name, ' part2']);
%     figfile=fullfile(save_res_folder,[strrep(crosssings(i).animal_name,' ','_'),'_linetraj_2']);
% %     saveas(fig,[figfile,'.fig']);
%     saveas(fig2,[figfile,'.png']);
%     saveas(fig2,[figfile,'.svg']);
%     saveas(fig2,[figfile,'.pdf']);
% end