function analyze_random_vertical_bar_opto_recording_subsmaple(red_signal, opto_burst_cont, sample_rate, logfile,resroot, depth_tip)
%reconstructs the stimulus of the random vertical bars, finds the time
%stamps of each frame and correlates the spiking activity with the stimulus

probe_length=800; %um
if ~exist('depth_tip','var') || isempty(depth_tip)
    depth_tip=probe_length;
end

[folder_to_analyze, log_name,~]  = fileparts(logfile);

%% reconstruct the bar locations
stimparam=load(logfile);
[bar_loc,barwidth_pix]  = vertical_random_bar_optointervals_reconstruct(logfile);


%% clean up red frames
sec_per_frame=0.0167;
frame_duration = int32(sec_per_frame*sample_rate);

redsum = movsum(red_signal, [0 frame_duration]);
rf=redsum>10;
red_signal_clean = rf & red_signal; %this will remove spurious red signal in between red frames

%time of red frames
rft=find(red_signal_clean);

%remove red frames registered within <2 frames
no_doubles = floor(2*frame_duration);
for i=1:length(rft)
    red_signal_clean(rft(i)+1:rft(i)+no_doubles)=0;
end

%clean red frames and get time stamp of each frame
[frame_timing, red_signal_i] = get_clean_frame_timing(red_signal_clean, logfile);

%% create an array with frame ids
% for each sample point generate the index of the frame presented
frame_index = zeros(size(red_signal_i));
nredframes=sum(red_signal_i);
framinds=1:5:5*nredframes;
redframesi=find(red_signal_i);
%first red frame
ist=redframesi(1);
%last red frame
ien=redframesi(end);
frame_index(ist:ien)=interp1(redframesi,framinds,ist:ien,'linear','extrap');
%round down values to get frame indices
frame_index=floor(frame_index);

%estimate length of the stimuli with opto
lenstim_opto_samples=sum(opto_burst_cont);
lenstim_opto_secs=lenstim_opto_samples/sample_rate;

%remove 0.1sec after opto to remove burst of neurons on opto release
optooff_start=find(diff([opto_burst_cont,0])==-1);
patw=round(0.1*sample_rate);
paddarray=zeros(size(opto_burst_cont));
mxN=numel(paddarray);
for i=1:length(optooff_start)
    eni=min(mxN,optooff_start(i)+patw);
    paddarray(optooff_start(i):eni)=1;
end
%estimate length of the stimuli without opto
lenstim_noopto_samples=ien-ist+1-lenstim_opto_samples - sum(paddarray);
lenstim_noopto_secs=lenstim_noopto_samples/sample_rate;

%collect all spike info files automatically
spike_cluster_file= fullfile(folder_to_analyze,'spike_clusters.npy');
spike_time_file= fullfile(folder_to_analyze,'spike_times.npy');
cluster_info_file= fullfile(folder_to_analyze,'cluster_group.tsv');
cluster_info_gen_file= fullfile(folder_to_analyze,'cluster_info.tsv');

[~,folder_name,~]=fileparts(folder_to_analyze);
res_folder= fullfile(resroot, folder_name);
rffigname=fullfile(resroot,[folder_name,'_rffigs.pdf']);
rfresfile=fullfile(resroot,[folder_name,'_rfinfo.mat']);


%     %make the folder with results
%     if exist(res_folder,'dir')==0
%         mkdir(res_folder);
%     end
%
%get depth info
ids_depth = get_clusters_depth(cluster_info_gen_file);

%frame_history_window
frame_history=30;
% frames into future
frames_future=5;

%cluster label and time of each spike
cluster_id = readNPY(spike_cluster_file);
spike_t = readNPY(spike_time_file);

%find ids of good clusters
goodclusters = get_clusters_by_property(cluster_info_file,'good');
muaclusters = get_clusters_by_property(cluster_info_file,'mua');
goodclusters=[goodclusters,muaclusters];

ncl=size(ids_depth,1);
RF(ncl) = struct();
rfstack=zeros(frame_history+frames_future, stimparam.tex_width);
sematrix=ones(1,barwidth_pix);
RFi=zeros(frame_history+frames_future,stimparam.tex_width);
%% for each neuron compute RF with/without opto
for ci=1:ncl
    %cluster id as after kilosort
    idc=ids_depth(ci,1);
    RF(ci).id = idc;
    RF(ci).depth = depth_tip - probe_length + ids_depth(ci,2);
    if ismember(idc,muaclusters)
        RF(ci).property= 'mua';
    elseif ismember(idc,goodclusters)
        RF(ci).property= 'good';
    else
        RF(ci).property= 'noise';
    end

    %select spikes per cluster
    spikes=find(cluster_id==idc);

    %get the time of spikes
    ts=spike_t(spikes);

    %convert spikes into 0,1 array
    spike_temp_array=zeros(1,length(red_signal_i),'logical');
    spike_temp_array(ts)=1;

    %find spikes during opto
    spike_opto=spike_temp_array&opto_burst_cont;
    tso=find(spike_opto);
    numspikes_opto = length(tso);

    %find spikes without opto
    spike_noopto=spike_temp_array&~(opto_burst_cont|paddarray);
    tsno=find(spike_noopto);
    numspikes_noopto = length(tsno);

    numspikes=min(numspikes_opto, numspikes_noopto);
    %subsample spikes
    if numspikes_opto > numspikes
        idx=randperm(numspikes_opto, numspikes);
        tso= tso(idx);
    end
    if numspikes_noopto > numspikes
        idx=randperm(numspikes_noopto, numspikes);
        tsno= tsno(idx);
    end


    %% generate RF with opto
    RFi(:)=0;
    for is=1:numspikes
        frameid_at_spike=frame_index(tso(is));
        if frameid_at_spike<frame_history, continue, end
        if frameid_at_spike>length(bar_loc)-frames_future, continue, end
        frameid_start=frameid_at_spike-frame_history+1;
        inds=frameid_start:frameid_at_spike+frames_future;
        rfstack(:)=0;
        lininds=sub2ind(size(rfstack),1:frame_history+frames_future,bar_loc(inds));
        rfstack(lininds)=1;
        rfstack=imdilate(rfstack,sematrix);
        RFi=RFi+rfstack;
    end
    RF(ci).RFopto=RFi./numspikes;
    RF(ci).numspikes_opto=numspikes_opto;
    RF(ci).numspikes_opto_used=numspikes;
    RF(ci).fr_opto=numspikes_opto/lenstim_opto_secs;

    %% generate RF without opto
    RFi(:)=0;
    for is=1:numspikes
        frameid_at_spike=frame_index(tsno(is));
        if frameid_at_spike<frame_history, continue, end
        if frameid_at_spike>length(bar_loc)-frames_future, continue, end
        frameid_start=frameid_at_spike-frame_history+1;
        inds=frameid_start:frameid_at_spike+frames_future;
        rfstack(:)=0;
        lininds=sub2ind(size(rfstack),1:frame_history+frames_future,bar_loc(inds));
        rfstack(lininds)=1;
        rfstack=imdilate(rfstack,sematrix);
        RFi=RFi+rfstack;
    end
    RF(ci).RFnoopto=RFi./numspikes;
    RF(ci).numspikes_noopto=numspikes_noopto;
    RF(ci).numspikes_noopto_used=numspikes;
    RF(ci).fr_noopto=numspikes_noopto/lenstim_noopto_secs;

    %number of spikes during opto release
    RF(ci).numspikes_optorelease=sum(spike_temp_array&paddarray);
    RF(ci).fr_optorelease= RF(ci).numspikes_optorelease/(sum(paddarray)/sample_rate);
end

%% flip array to match mouse's view
%mouse looks towards the projector, hence left-right are flipped
for ci=1:ncl
    RFtemp=RF(ci).RFnoopto;
    RFtemp=fliplr(RFtemp);
    RF(ci).RFnoopto=RFtemp;

    RFtemp=RF(ci).RFopto;
    RFtemp=fliplr(RFtemp);
    RF(ci).RFopto=RFtemp;
end

%% find and crop the area were stimulus was presented
x1=Inf;
x2=1;
for ci=1:ncl
    maxval=max(RF(ci).RFnoopto(:));
    if isnan(maxval) || maxval==0
        continue;
    end
    varx=var(RF(ci).RFnoopto);
    x1=min(x1,find(varx>0,1,"first"));
    x2=max(x2,find(varx>0,1,"last"));
end
x1=x1+10;
x2=x2-10;

%% compute different noise estimates
for ci=1:ncl
    %skip clusters with no rf
    maxval=max(max(RF(ci).RFopto(:)),max(RF(ci).RFnoopto(:)));
    if isnan(maxval) || maxval==0
        RF(ci).snr_var_opto=nan;
        RF(ci).snr_var_noopto=nan;
        RF(ci).logp2n_opto=nan;
        RF(ci).logp2n_noopto=nan;
        continue;
    end
    RF(ci).snr_var_opto=snr_variance_RFestimation(RF(ci).RFopto(:,x1:x2), frame_history, 1);
    RF(ci).logp2n_opto=peak2noise_log(RF(ci).RFopto(:,x1:x2), frame_history);
    RF(ci).snr_var_noopto=snr_variance_RFestimation(RF(ci).RFnoopto(:,x1:x2), frame_history, 1);
    RF(ci).logp2n_noopto=peak2noise_log(RF(ci).RFnoopto(:,x1:x2), frame_history);
end


%% save images of the RFs
for ci=1:ncl
    if isnan(RF(ci).snr_var_opto)
        continue;
    end
    maxval=max(max(RF(ci).RFopto(:)),max(RF(ci).RFnoopto(:)));
    %     if maxval==0
    %         continue;
    %     end

    fig = figure;
    t=tiledlayout(1,2,'TileSpacing','compact');

    nexttile;
    imagesc(RF(ci).RFopto, [0,maxval]);
    colormap(parula);
    hold on;
    plot([0.5,size(RF(ci).RFopto,2)+0.5],[frame_history+0.5 frame_history+0.5],'--w');
    title({'RF with opto'; ['p2n=',num2str(RF(ci).logp2n_opto)];...
        ['snrvar=',num2str(RF(ci).snr_var_opto)]});

    nexttile;
    imagesc(RF(ci).RFnoopto, [0,maxval]);
    colormap(parula);
    hold on; plot([0.5,size(RF(ci).RFnoopto,2)+0.5],[frame_history+0.5 frame_history+0.5],'--w');
    title({'RF without opto'; ['p2n=',num2str(RF(ci).logp2n_noopto)];...
        ['snrvar=',num2str(RF(ci).snr_var_noopto)]});

    title(t,['Cluster ', num2str(ci), '; depth ',num2str(RF(ci).depth)]);
    exportgraphics(fig,rffigname, 'BackgroundColor','none','ContentType','vector','Append', true);
    close(fig);
end


%% save variables for the RF analysis
save(rfresfile, 'RF', 'x1', 'x2', 'frame_history', 'stimparam', '-v7.3');
end

