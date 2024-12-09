repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\TNT_Sham_Max';
repofolder2= '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\OpenField\OpenField5\DLC_cropped2_CSV';
save_res_folder = 'D:\data\tnt_data\traj';

fi=dir(fullfile(repofolder ,'*.csv'));
ki=1;
fl=dir(fullfile(repofolder ,'*.csv'));
reminds=[];
for i=1:length(fi)
    if startsWith(fi(i).name,'._')
        reminds=[reminds,i];
    end
end
fl(reminds)=[];
ki=length(fl)+1;
fi=dir(fullfile(repofolder2 ,'*.csv'));
for i=1:length(fi)
    if ~startsWith(fi(i).name,'._')
        fl(ki)=fi(i);
        ki=ki+1;
    end
end


% %
% % %padded platform border
% % bd=600;
% % flipped_platform=5:9;
% % platformbins=0:0.1:2;
% % platformbins_1s=0.1:0.1:1;
% platform_vals=[100:100:bd];
% nb=numel(platform_vals);
% % nbins1s=numel(platformbins_1s);

allpts={};
flipped_platform=[5:9];

for i=1:length(fl)
    filename=fullfile(fl(i).folder, fl(i).name);
    [~,animal,~]=fileparts(fl(i).name);
    animal_id=strrep(animal,'_',' ');

    opts = detectImportOptions(filename);
    if i<10
        opts.SelectedVariableNames = {'Var5','Var6','Var7','Var17','Var18','Var20','Var21','Var23','Var24','Var26','Var27'};
        posxy = readmatrix(filename, opts);
        posxy(1:2,:)=[];
        bottom_left_top_right=[median(posxy(:,9),'omitnan'), median(posxy(:,4),'omitnan'),...
            median(posxy(:,5),'omitnan'),median(posxy(:,6),'omitnan')];
        rear_xy=posxy(:,1:2);
        rear_xy_ll=posxy(:,3);
    else
         opts.SelectedVariableNames = {opts.VariableNames{5},opts.VariableNames{6}};
         rear_xy=readmatrix(filename, opts);
         if iscell(rear_xy)
          rear_xy = cellfun(@str2double, rear_xy);
         end
         xynan=isnan(rear_xy(:,1));
         nansum=movsum(~xynan,[0,100]);
         fn=find(nansum>100,1,"first");
         rear_xy(1:fn-1,:)=[];

         opts.SelectedVariableNames = opts.VariableNames{7};
         rear_xy_ll=readmatrix(filename, opts);
         rear_xy_ll(1:fn-1)=[];
    end

    %interpolate nan values in ear positions
    nfr=size(rear_xy,1);
    inds=1:nfr;
    naninds=isnan(rear_xy(:,1)) | rear_xy_ll<0.8;
    allvals=interp1(inds(~naninds),rear_xy(~naninds,1),inds,"linear","extrap");
    rear_xy(:,1)=allvals;   
    allvals=interp1(inds(~naninds),rear_xy(~naninds,2),inds,"linear","extrap");
    rear_xy(:,2)=allvals;

    if ismember(i,flipped_platform)
        maxx=max(rear_xy(:,1));
        rear_xy(:,1)=maxx-rear_xy(:,1);
    end

    indscol=[1:size(rear_xy,1)]/size(rear_xy,1);
    fig =figure; scatter(rear_xy(:,1),rear_xy(:,2),[],indscol,'.');
    title(animal_id(1:min(length(animal_id),50)));
    figfile=fullfile(save_res_folder,['allpts_',animal(1:min(length(animal),50))]);
    saveas(fig,[figfile,'.fig']);
    saveas(fig,[figfile,'.png']);
    saveas(fig,[figfile,'.svg']);

    allpts(i).animal_name=animal_id;
    allpts(i).xy=rear_xy;
end

save(fullfile(save_res_folder,'allpts.mat'),"allpts");


%     if ismember(i,flipped_platform)
%         maxx=max(rear_xy(:,1));
%         rear_xy(:,1)=maxx-rear_xy(:,1);
%     end

%     crosssings(i).animal_name=animal_id;
%     crosssings(i).bottom_left_top_right=bottom_left_top_right;
%     crosssings(i).xy=rear_xy;


    %
    % xlim(bottom_left_top_right([2,4]));
    % ylim(bottom_left_top_right([1,3]));
    %
    % title(strrep(['exits for ',fl(i).name],'_',' '));
    % figfile=fullfile(save_res_folder,['exits_',crosssings(i).animal_name]);
    % saveas(fig,[figfile,'.fig']);
    % saveas(fig,[figfile,'.png']);
    % saveas(fig,[figfile,'.svg']);
% end
