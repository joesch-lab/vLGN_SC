% repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\DLC_Cliff\TeLC_weakVLGN_ZI';
repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\DLC_Cliff\Field3_Crpp_DLC_analysis';
anikeep={'714', '717', '718', '721'};
save_res_folder = 'D:\data\resub\data\cliff\medial_thalamus';
fps=60;

% load(fullfile(save_res_folder ,'arena_info_updated.mat'));

addlist=dir(fullfile(repofolder,'*.csv'));

allpts={};
ilen=length(allpts);
maxTplot=10*60*fps; %min

ki=1;
for k=1:length(addlist)
    %find matching mp4 file

    if ~contains(addlist(k).name, anikeep)
        continue;
    end
    disp(addlist(k).name);

    currfile=addlist(k).name;
    startname=currfile(1:strfind(currfile,'cropped')-1);
    vidfile=dir(fullfile(repofolder, [startname, '*cropped.mp4']));
    vidfile_name=vidfile.name;
    animal=vidfile_name(1:end-length('cropped.mp4'));
    animal_id=strrep(animal,'_',' ');
    vidfile=fullfile(vidfile.folder, vidfile.name);
   
    

    filename=fullfile(repofolder, currfile);
   
    opts = detectImportOptions(filename);
    if any(contains(opts.VariableNames,'Var')) %i<10
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

        opts.SelectedVariableNames = opts.VariableNames{7};
        rear_xy_ll=readmatrix(filename, opts);
        if iscell(rear_xy_ll)
            rear_xy_ll = cellfun(@str2double, rear_xy_ll);
        end

        xynan=isnan(rear_xy(:,1)) | rear_xy_ll<0.9;
        nansum=movsum(~xynan,[0,100]);
        fn=find(nansum>100,1,"first");
        rear_xy(1:fn-1,:)=[];
        rear_xy_ll(1:fn-1)=[];
    end
   
    %show the image and pts, allow user to define how to crop
    vr= VideoReader(vidfile);    
    im=readFrame(vr);
    cm=jet(size(rear_xy,1));
    for ii=1:size(rear_xy,1)
        if any(isnan(rear_xy(ii,:))), continue; end;
        im(round(rear_xy(ii,2)),round(rear_xy(ii,1)),:)=uint8(cm(ii,:)*255);
    end
    
    figure('Name', 'Select corners: left_top, left_bottom, right_bottom, right_top, left_top');
    [BW,xi2,yi2] = roipoly(im);
    
    xmina=mean(xi2(1:2));
    xmaxa=mean(xi2(3:4));
    ymina=mean(yi2([1,4]));
    ymaxa=mean(yi2(2:3));

    %find points with low likelihood or outside the platform
    nfr=size(rear_xy,1);
    inds=1:nfr;
    naninds=isnan(rear_xy(:,1)) | rear_xy_ll<0.9;
    nani=[find(rear_xy(:,1)<xmina);find(rear_xy(:,1)>xmaxa);...
      find(rear_xy(:,2)<ymina);find(rear_xy(:,2)>ymaxa)];
    nani=unique(nani(:));
    naninds(nani)=1;
    rear_xy(naninds,:)=NaN;

    %interpolate nan values
    allvals=interp1(inds(~naninds),rear_xy(~naninds,1),inds,"nearest");
    rear_xy(:,1)=allvals;
    allvals=interp1(inds(~naninds),rear_xy(~naninds,2),inds,"nearest");
    rear_xy(:,2)=allvals;


    %rescale all values to [1..1000]
    xmin=min(rear_xy(:,1));
    ymin=min(rear_xy(:,2));
    xmax=max(rear_xy(:,1));
    ymax=max(rear_xy(:,2));
    rear_x_sc=(rear_xy(:,1)-xmin+1)/(xmax-xmin+1)*1000;
    rear_y_sc=(rear_xy(:,2)-ymin+1)/(ymax-ymin+1)*1000;

    rear_xy(:,1)=rear_x_sc;
    rear_xy(:,2)=rear_y_sc;

   maxTplot_i=min(maxTplot, size(rear_xy,1));
    indscol=[1:maxTplot_i]/maxTplot_i;

    fig =figure; scatter(rear_xy(1:maxTplot_i,1),rear_xy(1:maxTplot_i,2),[],indscol,'.');
    title(animal_id(1:min(length(animal_id),50)));
    figfile=fullfile(save_res_folder,['allpts_',animal(1:min(length(animal),50))]);
    saveas(fig,[figfile,'.fig']);
    saveas(fig,[figfile,'.png']);
    saveas(fig,[figfile,'.svg']);

    fig =figure; plot(rear_xy(1:maxTplot_i,1),rear_xy(1:maxTplot_i,2),'k');
    title(animal_id(1:min(length(animal_id),50)));
    figfile=fullfile(save_res_folder,['allpts_lin_',animal(1:min(length(animal),50))]);
    saveas(fig,[figfile,'.fig']);
    saveas(fig,[figfile,'.png']);
    saveas(fig,[figfile,'.svg']);

    allpts(ilen+ki).animal_name=animal_id;
    allpts(ilen+ki).frame_start=fn;
    allpts(ilen+ki).xy=rear_xy;

     arena_info(ilen+ki).corners=[xi2(1), xi2(2), xi2(3), xi2(4), yi2(1), yi2(2), yi2(3), yi2(4)];    
     arena_info(ilen+ki).dlc_name=currfile;
     arena_info(ilen+ki).vid_name=vidfile_name;
     arena_info(ilen+ki).platform_left=1;
     ki=ki+1;
end

save(fullfile(save_res_folder,'new_allpts_updated.mat'),"allpts");
save(fullfile(save_res_folder,'arena_info_updated.mat'),"arena_info");

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
