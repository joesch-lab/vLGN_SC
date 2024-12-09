% loads csv file with the location of body parts annotated by DeepLabCut
% crops the values within the arean, plots traces of each recording

% repofolder='\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\OpenField\all_cropped_vid';
% repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\OpenField\all_dlc_files';
% repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\TNT_Sham_Max';
% repofolder2= '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\OpenField\OpenField5\DLC_cropped2_CSV';

repofolder = '\\fs.ista.ac.at\dfsgroup\joeschgrp\VZ\OpenField\Sham_NewJan2024\DLC';
save_res_folder = 'D:\data\tnt_data\traj\all';

fps=60;
plot_traces=false; %set true for ploting traces of recordings

load(fullfile(save_res_folder ,'arena_info_updated.mat'));
addlist=dir(fullfile(repofolder,'*.csv'));
if exist(fullfile(save_res_folder,'new_allpts.mat'),"file")~=0
    load(fullfile(save_res_folder,'new_allpts.mat'));
else
    allpts={};
end

ilen=length(allpts);
%% iterate trhu all recordings
for k=1:length(addlist)
    %find matching mp4 file
    currfile=addlist(k).name;
    startname=currfile(1:strfind(currfile,'cropped')-1);
    vidfile=dir(fullfile(repofolder, [startname, '*cropped.mp4']));
    vidfile_name=vidfile.name;
    animal=vidfile_name(1:end-length('cropped.mp4'));
    animal_id=strrep(animal,'_',' ');
    vidfile=fullfile(vidfile.folder, vidfile.name);

    %import dlc labels
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

        xynan=isnan(rear_xy(:,1)) | rear_xy_ll<0.9; %find non-existent or low-likelihood labels
        nansum=movsum(~xynan,[0,100]);
        fn=find(nansum>100,1,"first"); %find first high-likelihood labels: start of the experiment
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
    figure, [BW,xi2,yi2] = roipoly(im);
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

    % plot traces
    if plot_traces
        t15=min(15*60*fps, size(rear_xy,1));
        indscol=[1:t15]/t15;

        fig =figure; scatter(rear_xy(1:t15,1),rear_xy(1:t15,2),[],indscol,'.');
        title(animal_id(1:min(length(animal_id),50)));
        figfile=fullfile(save_res_folder,['allpts_',animal(1:min(length(animal),50))]);
        saveas(fig,[figfile,'.fig']);
        saveas(fig,[figfile,'.png']);
        saveas(fig,[figfile,'.svg']);

        fig =figure; plot(rear_xy(1:t15,1),rear_xy(1:t15,2),'k');
        title(animal_id(1:min(length(animal_id),50)));
        figfile=fullfile(save_res_folder,['allpts_lin_',animal(1:min(length(animal),50))]);
        saveas(fig,[figfile,'.fig']);
        saveas(fig,[figfile,'.png']);
        saveas(fig,[figfile,'.svg']);
    end

    %save

    allpts(ilen+k).animal_name=animal_id;
    allpts(ilen+k).frame_start=fn;
    allpts(ilen+k).xy=rear_xy;

    arena_info(ilen+k).corners=[xi2(1), xi2(2), xi2(3), xi2(4), yi2(1), yi2(2), yi2(3), yi2(4)];
    arena_info(ilen+k).dlc_name=currfile;
    arena_info(ilen+k).vid_name=vidfile_name;
    arena_info(ilen+k).platform_left=1;
end

save(fullfile(save_res_folder,'new_allpts_updated.mat'),"allpts");
save(fullfile(save_res_folder,'arena_info_updated.mat'),"arena_info");