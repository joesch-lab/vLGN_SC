% script to collect params of all random bar stimulus: width of the bar and
% opto burst cycle duration
clear all;

folderall='D:\data\VerticalBars_wOpto';
dirlist=dir(folderall);
barwidt=[];
opto_on_off=[];
for di=1:length(dirlist)
    if dirlist(di).isdir==0
        continue;
    end
    subdirname=fullfile(folderall,dirlist(di).name);
    paramfile=dir(fullfile(subdirname,'stimuli*log.mat'));
    if isempty(paramfile)
        continue;
    end
    paramfile=fullfile(paramfile.folder,paramfile.name);
    params=load(paramfile);
    barwidt=[barwidt,params.barwidth_deg];
    opto_on_off=[opto_on_off;[params.period_on, params.period_off]];
end