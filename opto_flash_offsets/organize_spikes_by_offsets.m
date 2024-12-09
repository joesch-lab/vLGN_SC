function organize_spikes_by_offsets(data, plot_histogram_yn, save_folder)

% save_folder = 'D:\NewMap\offsets_res';
[~,exp_folder,~]=fileparts(data.exp_info.folder);
im_save_folder =fullfile(save_folder,'cluster_im');
if ~exist(im_save_folder,'dir')
    mkdir(im_save_folder);
end

SR=data.trigger_data.sampling_rate;
N=length(data.trigger_data.red_signal);

%get the on edge of each visual flash
[on_vis,~] = get_signal_on_off_edge(data.trigger_data.red_signal);
[on_opto,~] = get_signal_on_off_edge(data.trigger_data.opto_cont);

%get params of flash and opto
%duration of each bout (no other flashes within this interval)
bout_t  = round(data.stimparams.bout_len*SR);
bout_t2=round(bout_t/2);

%duration of each flash
tflash_samples  = round(data.stimparams.flash_t*SR);
%duration of opto
topto_samples  = round(data.stimparams.opto_dur*SR);
%all offset types
offsets = data.stimparams.opto_offsets;

%for each flash get start and end of the bout and the offset type
nreps=length(on_vis);
n_offs = length(offsets);
ntrials=round(nreps/n_offs);
bout_st_en=zeros(nreps,2);
offset_sec=zeros(nreps,1);
flash_offset_trial = zeros(nreps,3);
trial_counter=zeros(n_offs ,1);

trial_array=zeros(1,N);
for i=1:nreps
    bout_st_en(i,1)=on_vis(i)-bout_t2;
    bout_st_en(i,2)=on_vis(i)+bout_t2;
    %find opto ON edge within [bout_st_en(i,1), bout_st_en(i,2)]
    opto_i=find((on_opto>=bout_st_en(i,1))&(on_opto<=bout_st_en(i,2)));
    opto_t=on_opto(opto_i);
    dt=(opto_t-on_vis(i))/SR;
    [~,io] = min(abs(offsets-dt));
    offset_sec(i)=offsets(io);
    trial_counter(io)=trial_counter(io)+1;
    flash_offset_trial(i,:)=[i,io,trial_counter(io)];
    trial_array(bout_st_en(i,1):bout_st_en(i,2))=i;
end

%for each cluster organize spikes into the flash bouts
ncl=length(data.spike_data);
%init data structure

offsets_trials=zeros(1,n_offs);
for ofi =1:n_offs
    offsets_trials(ofi)=max(flash_offset_trial(flash_offset_trial(:,2)==ofi,3));
end

flash_spikes={};
all_depth_CS=[];
for ci=1:ncl
    for ofi =1:n_offs
        for repi=1:offsets_trials(ofi)
            flash_spikes(ci).offset(ofi).rep(repi).sp=[];
        end
    end
    all_depth_CS=[all_depth_CS,data.spike_data(ci).SC_depth];
end

for ci=1:ncl
    tsp = single(data.spike_data(ci).spike_timing);
    for ti=1:length(tsp)
        flash_id = trial_array(tsp(ti));
        if flash_id==0, continue; end
        offset_i=flash_offset_trial(flash_id,2);
        rep_i=flash_offset_trial(flash_id,3);
        start_bout_i=bout_st_en(flash_id,1);
        flash_spikes(ci).offset(offset_i).rep(rep_i).sp=[flash_spikes(ci).offset(offset_i).rep(rep_i).sp,tsp(ti)-start_bout_i+1];
    end
end

%make plot per cluster
for ci=1:ncl
    f=figure;
    for ofi =1:n_offs
        for repi=1:offsets_trials(ofi)
            row_i = sum(offsets_trials(1:ofi-1))+repi;
            ts=flash_spikes(ci).offset(ofi).rep(repi).sp/SR;
            yvals=repelem(row_i,length(ts));
            plot(ts,yvals,'.k');
            hold on;
        end
    end
    for ofi =1:n_offs
        hold on;
        st_t=(bout_t2+offsets(ofi)*SR)/SR;
        en_t=st_t+topto_samples/SR;
        r1=sum(offsets_trials(1:ofi-1))+0.5;
        r2=sum(offsets_trials(1:ofi))+0.5;
        p1 = patch([st_t,en_t,en_t,st_t],[r1, r1, r2,r2],[0,0,1]);
        p1.EdgeColor = 'none';
        p1.FaceAlpha = 0.2;
        hold on;
    end
    st_t=bout_t2/SR;
    en_t=st_t+tflash_samples/SR;
    p2 = patch([st_t,en_t,en_t,st_t],[0.5, 0.5, sum(offsets_trials)+0.5,sum(offsets_trials)+0.5],[0,1,0]);
    p2.FaceAlpha = 0.2;
    p2.EdgeColor = 'none';
    yticks([]);
    ylim([0,sum(offsets_trials)+1]);
    xlim([0,bout_t/SR]);
    title({['cluster ',num2str(ci),'; depth ',num2str(data.spike_data(ci).SC_depth)]});
    filename = ['cluster_',num2str(ci)];
    % savefile=fullfile(save_folder,[filename,'.fig']);
    % savefig(f,savefile);
    savefile=fullfile(im_save_folder,[filename,'.pdf']);
    exportgraphics (f, savefile, 'BackgroundColor','none','ContentType','vector');
    savefile=fullfile(im_save_folder,[filename,'.png']);
    saveas(f,savefile);
    close(f);
end


%for each cluster, flash offset and trial make a FR histograms,
% average across trials
binlen=0.005; %5 ms
bine = 0:binlen:(bout_t/SR);
binc = bine(1:end-1)+binlen/2;
nbins=length(bine)-1;
hist_i=zeros(ncl,n_offs,nbins);
%responses to on edge of the flash
flashon_st_bin=floor((bout_t2/SR)/binlen);
flashon_en_bin=ceil((bout_t2/SR+0.1)/binlen);
flashoff_st_bin=floor(((bout_t2+tflash_samples)/SR)/binlen);
flashoff_en_bin=ceil(((bout_t2+tflash_samples)/SR+0.1)/binlen);

%start and end of opto in bins
opto_st_t=(bout_t2/SR+offsets)/binlen;
opto_en_t=opto_st_t+topto_samples/SR/binlen;

on_resp=zeros(ncl,n_offs);
off_resp=zeros(ncl,n_offs);
maxflash_resp=zeros(ncl,n_offs);
for ci=1:ncl
    for ofi =1:n_offs
        tsall=[];
        for repi=1:offsets_trials(ofi)
            tsall=[tsall,flash_spikes(ci).offset(ofi).rep(repi).sp/SR];            
        end
        hist_i(ci,ofi,:)=reshape(histcounts(tsall,bine),1,1,[]);
        hist_i(ci,ofi,:)=hist_i(ci,ofi,:)./offsets_trials(ofi); %average
    end
    hist_i(ci,:,:)=hist_i(ci,:,:)/binlen; %convert FR to spikes/sec

    on_resp(ci,:)=max(hist_i(ci,:,flashon_st_bin:flashon_en_bin),[],3);
    off_resp(ci,:)=max(hist_i(ci,:,flashoff_st_bin:flashoff_en_bin),[],3);
    maxflash_resp(ci,:)=max(hist_i(ci,:,flashon_st_bin:flashoff_en_bin),[],3);



    %make a figure per cluster
    if plot_histogram_yn
        f=figure;
        t=tiledlayout(3,2);
        nexttile([3,1]);
        imagesc(squeeze(hist_i(ci,:,:)));
        set(gca,'YDir','normal');
        colormap(flipud(colormap('gray')));
        hold on;

        st_t=(bout_t2/SR)/binlen;
        en_t=st_t+tflash_samples/SR/binlen;
        p = patch([st_t,en_t,en_t,st_t],[0.5, 0.5, n_offs+0.5,n_offs+0.5],[0,1,0]);
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.05;
        for ofi =1:n_offs
            hold on;
            p1 = patch([opto_st_t(ofi),opto_en_t(ofi),opto_en_t(ofi),opto_st_t(ofi)],[ofi-0.5, ofi-0.5, ofi+0.5,ofi+0.5],[0,0,1]);
            p1.EdgeColor = 'none';
            p1.FaceAlpha = 0.05;
        end
        colorbar();
        yticks([]);
        xticks(1:50:nbins);
        xticklabels(round(binc(1:50:end),1));


        nexttile;
        plot(offsets,on_resp(ci,:));
        ymax=1.2*max([max(on_resp(ci,:)),max(off_resp(ci,:)),max(maxflash_resp(ci,:))]);
        if ymax==0, ymax=1; end
        ylim([0,ymax]);
        ylabel('ON');
        nexttile;
        plot(offsets,off_resp(ci,:));
        ylim([0,ymax]);
        ylabel('OFF');
        nexttile;
        plot(offsets,maxflash_resp(ci,:));
        ylim([0,ymax]);
        ylabel('ON-OFF');

        title(t,{['cluster ',num2str(ci),'; depth ',num2str(data.spike_data(ci).SC_depth)]});
        filename = ['cluster_',num2str(ci),'_histFR'];
        savefile=fullfile(im_save_folder,[filename,'.fig']);
        savefig(f,savefile);
        savefile=fullfile(im_save_folder,[filename,'.png']);
        saveas(f,savefile);
        close(f);
    end
end


of_data.folder = exp_folder;
of_data.flashI = 255*(max(data.stimparams.flash_color)>10);
of_data.spikes = flash_spikes;
of_data.cluster_depth = all_depth_CS;
of_data.offsets = offsets;
of_data.boutlen_samples = bout_t;
of_data.flash_st = bout_t2;
of_data.flash_dur = tflash_samples;
of_data.opto_dur = topto_samples;
of_data.sr= SR;

of_data.histograms.hist = hist_i;
of_data.histograms.on_resp=on_resp;
of_data.histograms.off_resp=off_resp;
of_data.histograms.onoff_rest=maxflash_resp;
of_data.histograms.onoff_rest=maxflash_resp;
of_data.histograms.binlen=binlen;
of_data.histograms.binc=binc;
of_data.histograms.bin_flash_st=flashon_st_bin;
of_data.histograms.bin_flash_en=flashoff_st_bin;
of_data.histograms.bin_opto_st=opto_st_t;
of_data.histograms.bin_opto_en=opto_en_t;
savefile=fullfile(save_folder,'opto_flash_data.mat');
save(savefile,'of_data','-v7.3');
