% compute_flash_optooffsets_curve

%for all neurons from all recordings select neurons that are both opto and
%visual (after alreading running permutation and zeta tests)

%for the selected neurons compute the spike count and firing rates at different opto
%offsets, normalize the response by max firing rate, make significance
%tests and plot

ii=1;
curve_all=[];
spikes_before=[];
spikes_opto=[];

data_folder='D:\data\resub\data';
path_to_figures= 'D:\data\resub\figs';
resfolder= data_folder;

list_process='D:\data\resub\data\with_laser_list.txt';

% list_process='D:\NewMap\all_rec_list.txt';
folderlist = get_list_from_file(list_process);

nf=length(folderlist);

pthresh=0.01;

sr=20000;
r_win_s=0.2;
tbefore_s = r_win_s/2;
bin_dur_s=0.005;
bine=-tbefore_s:bin_dur_s:r_win_s;

ii=1;

stimdur=[]; %1*sr; % filter recordings by stim duration
nrec=0;
%% loop for partially pre-processed files
for fi=1:nf     
    fileid=folderlist{fi};        
    data_dir=fullfile(data_folder,fileid);
    ofdata_file = fullfile(data_dir,'opto_flash_data.mat');        
    load(ofdata_file);
    
    nrec=nrec+1;

    data_dir=fullfile(resfolder,fileid);
    param_file = fullfile(data_dir,'params_flash_opto_responses_all.mat');        
    load(param_file);

    ncl=length(params_vals);
    nofi=length(of_data.offsets);

    %make the mask to compute flash spikes during different offsets
    %duration of the bout
    tbout=of_data.boutlen_samples;
    %start, duration, end of the flash
    tf_st=of_data.flash_st;
    tf_dur = of_data.flash_dur;
    tf_en=tf_st+tf_dur-1;      
       
    offsets_hist=zeros(nofi,length(bine)-1);    
    spikes_200ms=zeros(1,nofi);
    for ci=1:ncl
        if ~(is_visual(params_vals(ci),pthresh)&& is_opto(params_vals(ci),pthresh))
            continue;
        end

        %get spikes during the flashes 
        offsets_hist(:)=0;
        spikes_200ms(:)=0;               
        ri=1;
        for ofi=1:nofi        
            numspi=0;
            spikes=[];
            nrep_i=size(of_data.spikes(ci).offset(ofi).rep,2);        
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp; 
                spi=(spi-tf_st)/sr;
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                spikes=[spikes,spi]; 
                numspi=numspi+numel(spi(spi>=0));
            end

            spikes_200ms(ofi)=numspi/nrep_i;
            offsets_hist(ofi,:)=histcounts(spikes,bine); 
        end        

        curve_all(ii).name=fileid;
        curve_all(ii).flash_dur = of_data.flash_dur;
        curve_all(ii).flash_I = of_data.flashI;
        curve_all(ii).cl_num=ci;
        curve_all(ii).depth= of_data.cluster_depth(ci);
        curve_all(ii).flash_hist=offsets_hist;         
        curve_all(ii).spike_count_200ms = spikes_200ms;
        ii=ii+1;
    end
end

disp(['number of recordings: ',num2str(nrec)]);


%% flash opto curve:
% for each neuron and offset plot responses

%select recordings of the same length
duration_selected=0.2*sr;%1*sr;%
flashI_selected=[]; %% all intensities

%% plot for sSC
make_flash_opto_curve_plot_spike_count(curve_all, duration_selected, flashI_selected, 0,500,'sSC',path_to_figures,sr);
%% plot for dSC
make_flash_opto_curve_plot_spike_count(curve_all, duration_selected, flashI_selected, 500,Inf,'dSC',path_to_figures,sr);
%% plot for all
make_flash_opto_curve_plot_spike_count(curve_all, duration_selected, flashI_selected, 0,Inf,'all',path_to_figures,sr);

% save raw data
save(fullfile(path_to_figures,'flash_opto_curve_spikecount_all_data.mat'),'curve_all','bine', 'pthresh');


function make_flash_opto_curve_plot_spike_count(curve_all, duration_selected, flashI_selected, minD,maxD,labelstr,path_to_figures,sr)
% flash opto curve:
% for each neuron and offset plot responses
% preselect neurons from recording of the specific length duration_selected
% preselect neurons from the given depth range [minD,maxD)

curvesall=[];
rec_name={};
for i=1:length(curve_all)
    if ~isempty(duration_selected) && curve_all(i).flash_dur~=duration_selected
        continue;
    end
    if ~isempty(flashI_selected) && curve_all(i).flash_I~=flashI_selected
        continue;
    end
    if ~(curve_all(i).depth >=minD && curve_all(i).depth<maxD)
        continue;
    end
    curvesall=cat(1,curvesall,curve_all(i).spike_count_200ms);
    rec_name=cat(1,rec_name,curve_all(i).name);
end

nall=size(curvesall,1);
rec_list=unique(rec_name);
nrec=numel(rec_list);
a_list=animal_list(rec_list);
nanimals=numel(a_list);
substr = ['n=', num2str(nall), ', nrec=', num2str(nrec), ', nanimals=', num2str(nanimals)];

maxvals=max(curvesall,[],2);
mean_curve=mean(curvesall./maxvals,'omitnan');

%% for each offset check the difference with the offset durtherst from opto
alpha=0.01;
nofi=size(curvesall,2);
curvesall_n=curvesall./maxvals;
flash_r=curvesall_n(:,end);
for oi=nofi-1:-1:1
    flash_ofi=curvesall_n(:,oi);

    [p,h,stats]  = signrank(flash_r,flash_ofi,'tail','both','alpha',alpha);
    signrank_testres_offsetcurve(oi).significance_level=alpha;
    signrank_testres_offsetcurve(oi).tail='both';
    signrank_testres_offsetcurve(oi).test_res_h=h;
    signrank_testres_offsetcurve(oi).pvalue=p;
    signrank_testres_offsetcurve(oi).zstatistics=stats.zval;
    signrank_testres_offsetcurve(oi).stats=stats;
    signrank_testres_offsetcurve(oi).spike_count_n=flash_ofi;
end
signrank_testres_offsetcurve(nofi).spike_count_n=flash_r;

%% figure
f=figure;
x = repelem(1:13, nall);
y = reshape(curvesall./max(curvesall,[],2),1,[]);
sch = swarmchart(x,y(:),15,'filled', 'MarkerFaceColor',[64/255, 149/255, 206/255], 'MarkerFaceAlpha', 0.5);
sch.XJitterWidth=0.7;
hold on;
plot(mean_curve,'Color',[64/255, 149/255, 206/255],'LineWidth',2)
ylim([0,1.2]);
% add asterics for significance of difference
for oi=nofi-1:-1:1
    hold on;
    if signrank_testres_offsetcurve(oi).pvalue<=0.001 % ***
        text(oi-0.25, 1.1, '***');
    elseif signrank_testres_offsetcurve(oi).pvalue<=0.01 % **
        text(oi-0.25, 1.1, '**');
    elseif signrank_testres_offsetcurve(oi).pvalue<=0.05 % *
        text(oi-0.25, 1.1, '*');
    else
        text(oi-0.25, 1.1, 'n.s.');
    end
end

xticks(1:2:13);
opto_offsets=-1.5:0.25:1.5;
xvals=fliplr(opto_offsets);
xticklabels(xvals(1:2:end)*duration_selected/sr);
xlabel('seconds before flash onset');

title({labelstr;['flash dur = ',num2str(duration_selected/sr),'s']; substr});
durms=round(duration_selected/sr*1000);
savefig(f,fullfile(path_to_figures,['flash_opto_curve_spikecount_',labelstr,'_dur',num2str(durms),'ms.fig']));
saveas(f,fullfile(path_to_figures,['flash_opto_curve_spikecount_',labelstr,'_dur',num2str(durms),'ms.svg']));
save(fullfile(path_to_figures,['flash_opto_curve_spikecount_',labelstr,'_dur',num2str(durms),'ms_data.mat']),...
    'curvesall_n','mean_curve', 'signrank_testres_offsetcurve','maxvals','duration_selected','opto_offsets',...
    'rec_list', 'nall', 'nrec', 'a_list', 'nanimals');
end

