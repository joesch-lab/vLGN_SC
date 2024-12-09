%% script analyzes neuronal responses to the black or white flashes. 
% The stimulus is 1 sec long flashes of various intencities (usually white,
% gray and black) randomly interleaved. The script extracts white
% (from black to white) or black (from white to black) flashes. PTHS, mean
% population responses and average firing rates per neuron are plotted

% O.Symonova
% 2024

% folder where all stimuli are located
repof='D:\data\tnt_data';
 
%% collect all files with linear step stimulus and extract flash responses
%colldect all '_data.mat' files
%check if they have linear flash stimulus
%collect recordings with neuronal activity
stimname='stimuli_linear_intensity';
datafilelist = collect_processed_stimuli_files(repof,stimname, 0, 1, 0);

binw=0.01; %bin width in secs
tbefore=0.5; %time to consider before the flash onset
tafter=1; %time after the flash onset
idblack=1; %number of the black flash in the stimulus file
idwhite=3; %number of the white flash in the stimulus file

%select to analyze black or white flashes
blackflash=0;
if blackflash
    resfolder = 'D:\data\resub\data\linear_step_flash\black_flash_short';       
    flash1=idwhite; %start from white  
    flash2=idblack; %black flash
else %white flash
    resfolder = 'D:\data\resub\data\linear_step_flash\white_flash_short';
    flash1=idblack; %start from black 
    flash2=idwhite; %white flash
end


linflash_datafile=fullfile(resfolder,'linflash_neuronalresponse.mat');
params_da={};
params_da.tbefore_s=tbefore;
params_da.tafter_s=tafter;
params_da.binw=binw;
params_da.idblack=idblack;
params_da.idwhite=idwhite;
params_da.flash1=flash1;
params_da.flash2=flash2;

k=1;
linstep_neurodata={};
for fi=1:length(datafilelist)
    clear data;
    datafile=datafilelist{fi};
    load(datafile);

    if isfield(data.exp_info, 'folder')
        ofi=data.exp_info.folder;
    else
        ofi=data.exp_info.fullpath;
    end
    aniid = get_animal_id_from_recordingname(ofi);
    if contains(datafile,'control',"IgnoreCase",true),  continue; end;

    [~,exp_name,~]= fileparts(datafile);
    if contains(exp_name,'tnt','IgnoreCase',true) %is it sham or telc animal
        istnt=true;
    else
        istnt=false;
    end

    %find id of the logfile with the stimulus
    logid = find(contains({data.exp_info.logfiles(:).name},stimname));
    if length(logid)>1
        disp(ofi);
    end
    logid=logid(1);
    %get the date of the experiment
    if isfield(data.exp_info.logfiles(logid), 'date')
        currdate=data.exp_info.logfiles(logid).date;
    else
        currdate=data.exp_info.logfiles(logid).file_tstamp;
    end
    stiparamfile = fullfile(ofi, data.exp_info.logfiles(logid).name);
    %load stimulus paremeters
    stimdata = load(stiparamfile);

    sr=data.trigger_data.sampling_rate;
    %start and end of the stimulus within the recording
    stimstart=data.exp_info.logfiles(logid).start_s;
    stimend=data.exp_info.logfiles(logid).end_s;

    spikes=data.spike_data;    

    %get the red signal during the recording and when each flash starts
    st=stimstart*sr;
    en=stimend*sr;
    stmin=max(1,round(st-params_da.tbefore_s*sr));
    enmax=min(length(data.trigger_data.red_signal), round(en+0.5*sr));
    rs=data.trigger_data.red_signal(stmin:enmax);
    %start of each new flash
    vis_on = get_vis_on_edge_no_other(rs,zeros(size(rs)), 0.25, 0.25, sr);
    vis_on = (vis_on+double(st-params_da.tbefore_s*sr)-1)/sr;

    %there are 6 types of flashes: 0->1; 1->0; 0.5->1; 1->0.5; 0->0.5; 0.5->0;
    %from a random sequence of pairs of flashes create one vector:  0->1->0.5->1->0.5... 
    flash_id=nan(1,numel(stimdata.present_order_all)*2);
    nreps=size(stimdata.present_order_all,1);
    nlevels=size(stimdata.present_order_all,2);
    for i=1:nreps
        for j=1:nlevels
            ij=(i-1)*nlevels*2+(j-1)*2;
            pair=stimdata.present_order_all(i,j);
            flash_id(ij+1)=stimdata.pairsI(pair,1);
            flash_id(ij+2)=stimdata.pairsI(pair,2);
        end
    end

    nflash=length(flash_id);
    %find all occurances of flashes [0,255,255]->[0,0,0]; id 3->1 or 3->1
   
    %for each flash type find reps
    flash_pairs=stimdata.pairsI;

    idp =[]; %indices of flashes

    id1=find(flash_id==flash1); %find id of the starting intensity
    id1=id1(id1<nflash);
    for j=1:length(id1)
        if flash_id(id1(j)+1)==flash2 %check if the nex flash matching the needed intensity
            idp = [idp, id1(j)];
        end
    end

    flash_st=vis_on(idp+1); %start time of each flash

    %get FR histograms on the flach
    [spike_hist, histcen] = get_spikes_on_times(flash_st, spikes, sr, tbefore, tafter, binw);

    %make zeta-test for the responsiveness to flashes (expect response within 0.4sec)
    neuron_zetavals = stim_resp_zetatest(flash_st,spikes, sr, 0.4);

    %depth of neurons from CSDA
    depthSC = [spikes(:).SC_depth]';

    %vector with the recording, animal id, istnt
    nn=size(spike_hist,1);
    recvec=ones(nn,1)*fi;
    aniid_fi=ones(nn,1)*aniid;
    istntvec=ones(nn,1)*istnt;

    linstep_neurodata(k).datafile=datafile;
    linstep_neurodata(k).aniid=aniid;
    linstep_neurodata(k).date=currdate;
    linstep_neurodata(k).istnt=istnt;
    linstep_neurodata(k).flash_colors=[stimdata.all_colors(flash1,:); stimdata.all_colors(flash2,:)];
    linstep_neurodata(k).stimstart=data.exp_info.logfiles(logid).start_s;
    linstep_neurodata(k).stimsend=data.exp_info.logfiles(logid).end_s;
    linstep_neurodata(k).flash_start=flash_st;
    linstep_neurodata(k).fr_hist = spike_hist;
    linstep_neurodata(k).fr_histcen = histcen;
    linstep_neurodata(k).neuron_zeta_pvals = [neuron_zetavals(:).z_pval]';
    linstep_neurodata(k).neuron_depth= depthSC;
    k=k+1;
end
save(linflash_datafile,'linstep_neurodata', 'params_da', '-v7.3');

%% transform structure into arrays for plotting
load(linflash_datafile);

histall=[];
tntvec=[];
depthvec=[];
zpvalvec=[];
anivec=[];
recvec=[];
for i=1:length(linstep_neurodata)
    histall=[histall;linstep_neurodata(i).fr_hist];
    depthvec=[depthvec;linstep_neurodata(i).neuron_depth];
    zpvalvec=[zpvalvec;linstep_neurodata(i).neuron_zeta_pvals];
    nn=size(linstep_neurodata(i).fr_hist,1);
    tntvec=[tntvec; ones(nn,1)*linstep_neurodata(i).istnt];
    anivec=[anivec; ones(nn,1)*linstep_neurodata(i).aniid];
    recvec=[recvec; ones(nn,1)*i];
end
histcen=linstep_neurodata(1).fr_histcen;


%% make figure of neuronal responses to flashes: PST histogram
minbinval=-0.1;
maxbinval=0.3;
minbin=find(histcen>=minbinval,1,"first");
maxbin=find(histcen<=maxbinval,1,"last");
histall=histall(:,minbin:maxbin);
histcen=histcen(:,minbin:maxbin);

idx=tntvec>0; %indices of TeLC animals
if blackflash
    titlestr={'TNT responses to black flash'};
    figsavename=fullfile(resfolder, 'neuronal_responses_flashes_black.pdf');
else
    titlestr={'TNT responses to white flash'};
    figsavename=fullfile(resfolder, 'neuronal_responses_flashes_white.pdf');
end

make_response_hist_figure(histall(idx,:), histcen, zpvalvec(idx), depthvec(idx), 0.01, -Inf, 500, titlestr, figsavename, 1, true);

idx=tntvec<1; %indices of sham animals
if blackflash, titlestr={'Ctrl responses to black flash'};
else, titlestr={'Ctrl responses to white flash'};
end
make_response_hist_figure(histall(idx,:), histcen, zpvalvec(idx), depthvec(idx), 0.01, -Inf, 500, titlestr, figsavename, 1, true);
save(fullfile(resfolder, 'neuronal_responses_allneurons.mat'), 'histall', 'histcen', 'zpvalvec', 'depthvec', 'tntvec', 'anivec', 'recvec');

%% mean population responses
idx1=tntvec<1;
idx2=tntvec>0;
if blackflash, titlestr='mean responses to OFF flashes';
else, titlestr='mean responses to ON flashes';
end
mean_population_plots(histall(idx1,:), histall(idx2,:), depthvec(idx1), depthvec(idx2), ...
    zpvalvec(idx1), zpvalvec(idx2),  anivec(idx1), anivec(idx2), recvec(idx1), recvec(idx2), ...
    histcen, -Inf, 500, 0.01, 0, titlestr, figsavename);

%% spike count and width fo responses
 quantify_response_width_spikecount(histall(idx1,:), histall(idx2,:), depthvec(idx1), depthvec(idx2), ...
    zpvalvec(idx1), zpvalvec(idx2),  anivec(idx1), anivec(idx2), recvec(idx1), recvec(idx2), ...
    histcen, -Inf, 500, 0.01, figsavename);
