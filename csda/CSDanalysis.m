function csda = CSDanalysis(data, csda_res_folder, make_figure_yn, vis_only, start_stim, end_stim, vis_on, dur_before_on_edge, dur_after_s, tbefore_mean)
% run the CSD analysis
% params:
% data: the structure with triggeres
% csda_res_folder: folder where to store results of the CSD
% make_figure_yn: 1 if make figures
% vis_only: 1 if run CSD only on visual flashes, 0 if also on opto flashes
%         start_stim/end_stim: start/end time in samples of the stimulus with flashes
%     vis_on: the time in samples of the start of the flashes if given
% dur_before_on_edge, dur_after_s: time in secs to consider before and after onset of the flash
% tbefore_mean: time in secs of the window to computer the mean before flash onset

if ~exist('vis_only',"var") || isempty(vis_only)
    vis_only=0;
end

if ~exist('start_stim',"var") || isempty(start_stim)
    start_stim=1;
end

if ~exist('end_stim',"var") || isempty(end_stim)
    end_stim=length(data.trigger_data.red_signal);
end

if ~exist('dur_after_s',"var") || isempty(dur_after_s)
    dur_after_s=0.2;
end

if ~exist('dur_before_on_edge',"var") || isempty(dur_before_on_edge)
    %max duration to avoid other signal depends on the stimulus
    % dur_before_on_edge= data.stimparams.bout_len/2;
    dur_before_on_edge= 0.5;
end

if ~exist('tbefore_mean',"var") || isempty(tbefore_mean)
    %max duration to avoid other signal depends on the stimulus
    % dur_before_on_edge= data.stimparams.bout_len/2;
    tbefore_mean= dur_before_on_edge;
end

redsiginal=data.trigger_data.red_signal(start_stim:end_stim);
sampling_rate = data.trigger_data.sampling_rate;
binfile=data.trigger_data.binfile;

%get the time of ON-OFF edge of viusal and opto flash
if  ~exist('vis_on',"var") || isempty(vis_on)    
    if ~vis_only
        opto_cont=data.trigger_data.opto_cont;
        vis_on   = get_vis_on_edge_no_other(redsiginal,opto_cont, dur_before_on_edge, dur_after_s, sampling_rate);
        opto_on  = get_opto_on_edge_no_other(redsiginal,opto_cont, dur_before_on_edge, dur_after_s, sampling_rate);
    else
        vis_on   = get_vis_on_edge_no_other(redsiginal,zeros(size(redsiginal)), dur_before_on_edge, dur_after_s, sampling_rate);
    end
end

%padding flashes: time before and after the start of the flash
tbefore = dur_before_on_edge*sampling_rate;
tafter=dur_after_s*sampling_rate;
nlen=tbefore+tafter+1;

%parameters of the silicone probe
% 32 channels and 25um between the channels
nchannels=32;
spacing=2.5e-5;

%read raw data from the silicone probe (before spike sorting)
fileid = fopen(binfile);

data_raw=int16(fread(fileid,'int16'));
fclose(fileid);
data_raw=reshape(data_raw,nchannels,[]);
data_raw=data_raw(:,start_stim:end_stim);
load('chanMapCamNeu_Lin-A32-25_August2022_5.mat');
data_raw=data_raw(flipud(chanMap)',:);
cha_labels=flipud(chanMap);
cha_labels=cha_labels-1;

%for the recordings on 220223 and 220315 the three deepest channels were
%corrupted, set the values to nan;
if isfield(data.exp_info, 'folder') && (contains(data.exp_info.folder,'220223') || contains(data.exp_info.folder,'220224') ||contains(data.exp_info.folder,'220315'))
    data_raw(end-2:end,:)=nan;
end

%remove high-freq signal by median filtering in 5ms window
win_noise=0.005*sampling_rate;

%for all ON visual edges run CSD analysis and take the mean
CSD_vis_on=zeros(nchannels,nlen);
enum=length(vis_on);
vi_von=0;
for i=1:enum
    istart=vis_on(i)-tbefore;
    iend = vis_on(i)+tafter;
    data2=single(data_raw(:,istart:iend)');
    data2= medfilt1(data2,win_noise,[],1);
    CSDoutput  = CSD_nopic(data2,sampling_rate,spacing);
    CSD_vis_on=CSD_vis_on+CSDoutput';
    vi_von=vi_von+1;
end
CSD_vis_on_av=CSD_vis_on/vi_von;

%for all ON opto edges run CSD analysis and take the mean
if ~vis_only
    CSD_optoon=zeros(nchannels,nlen);
    enum=length(opto_on);
    vi_oon=0;
    for i=1:enum
        istart=opto_on(i)-tbefore;
        iend = opto_on(i)+tafter;
        data2=single(data_raw(:,istart:iend)');
        data2= medfilt1(data2,win_noise,[],1);
        CSDoutput  = CSD_nopic(data2,sampling_rate,spacing);
        CSD_optoon=CSD_optoon+CSDoutput';
        vi_oon=vi_oon+1;
    end
    CSD_optoon_av=CSD_optoon/vi_oon;
end

%corrected channels of the most recent recordings, set to NaN
if isfield(data.exp_info, 'folder') && (contains(data.exp_info.folder,'220223') || contains(data.exp_info.folder,'220224') ||contains(data.exp_info.folder,'220315'))
    CSD_vis_on_av(end-3:end,:)=nan;
    if ~vis_only, CSD_optoon_av(end-3:end,:)=nan; end
end

%substract the mean value of the CSDA before the onset of the stimulus
tbefore_remove=max(1,tbefore-tbefore_mean*sampling_rate);
tbefore=tbefore_mean*sampling_rate;

%vis
CSD_vis_on_av=CSD_vis_on_av(:,tbefore_remove:end);
meantemp=mean(CSD_vis_on_av(:,1:tbefore),2);
CSD_vis_on_av=CSD_vis_on_av-meantemp;
inflection_channel_id_vis=find_inflection_point(CSD_vis_on_av);
var_t_vis = CSDA_response_timing(CSD_vis_on_av, tbefore+1, tbefore_mean, dur_after_s, sampling_rate); % get dynamics of csd
env_vis = csd_variance(CSD_vis_on_av,tbefore+1); %envelope

%opto
if ~vis_only
    CSD_optoon_av=CSD_optoon_av(:,tbefore_remove:end);
    meantemp=mean(CSD_optoon_av(:,1:tbefore),2);
    CSD_optoon_av=CSD_optoon_av-meantemp;
    inflection_channel_id_opto=find_inflection_point(-CSD_optoon_av);
    var_t_opto = CSDA_response_timing(CSD_optoon_av, tbefore+1, tbefore_s, dur_after_s, sampling_rate); % get dynamics of csd
    env_opto = csd_variance(CSD_optoon_av,tbefore+1); %envelope
end

%save data
csda.inflection_cha_above = cha_labels(inflection_channel_id_vis);
csda.channel_labels = cha_labels;
csda.CSD_vis_on = CSD_vis_on_av;
csda.CSD_vis_on_nrep = vi_von;
csda.var_t_vis = var_t_vis;
csda.envelope_vis = env_vis;
csda.flash_onset_samples = tbefore+1;
csda.cha_spacing_um = spacing*1e+6;
csda.binfile = binfile;

if ~vis_only
    csda.inflection_cha_above_opto = cha_labels(inflection_channel_id_opto);
    csda.CSD_opto_on = CSD_optoon_av;
    csda.CSD_opto_on_nrep = vi_oon;
    csda.var_t_opto = var_t_opto;
    csda.envelope_opto = env_opto;
end

if make_figure_yn
    csda_figure_all(csda,csda_res_folder, vis_only);
end
end

function inflection_id=find_inflection_point(csd)
%find min vals across time for each channel
min_scd=min(csd,[],2);
%find sink
[~,sink_id]=min(min_scd);
%find the time when the think is stronges
[~,t_sink]=min(csd(sink_id,:));
%find the last positive chanel before the sink at the time t_sink
inflection_id = find(csd(1:sink_id,t_sink)>0,1,'last');
end