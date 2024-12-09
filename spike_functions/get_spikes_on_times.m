function [spike_hist, binc] = get_spikes_on_times(tpts,spike_data, sr, tbefore, tafter, binw0)
% get spike histogram per unit, all reps combined together
% params:
% tpts: time points around which to compute histograms
% spike_data: data structure with spikes and unit information
% sr: sampling rate of spikes (20K for silicone probe, 50K for neuropixels)
% tbefore: time in secs to consider before event onset
% tafter: time in secs  to consider after event onset, for histogram computation
% binw0: width of the bin
%
% O.Symonova, 2024

nch=length(spike_data);
nt=length(tpts);
spikes_raw={};
for ci=1:nch
    if abs(sr-20000)<1
        spikest=double(spike_data(ci).spike_timing)/sr;
    else
        spikest=double(spike_data(ci).spike_timing);
    end
    for ti=1:nt
        tst=tpts(ti)-tbefore;
        ten=tpts(ti)+tafter;
        sp_sel=(spikest(spikest>=tst & spikest<=ten))-tpts(ti);
        spikes_raw(ci).rep(ti).spikes=sp_sel;
    end
end

%     for ci=1:nch
%         nsp=0;
%         for ti=1:nt
%             nsp=nsp+numel(spikes_raw(ci).rep(ti).spikes);
%         end
%         if nsp==0
%             continue;
%         end
%         f=figure;
%         for ti=1:nt
%             yval=ones(1,numel(spikes_raw(ci).rep(ti).spikes))*ti;
%             plot(spikes_raw(ci).rep(ti).spikes, yval, '.k');
%             hold on;
%         end
%         xlim([-tbefore,tafter]);
%         title(['cluster ',num2str(ci)]);
%     end
%
if ~exist('binw0','var') || isempty(binw0)
    nbin=20;
    binw=(tafter+tbefore)/nbin;
else
    binw=binw0;
    nbin=(tafter+tbefore)/binw;
end

bine=linspace(-tbefore,tafter,nbin+1);
binc=bine(1:end-1) + binw/2;

spike_hist=zeros(nch,nbin);
for ci=1:nch
    spi = [];
    for ti=1:nt
        spi=[spi,spikes_raw(ci).rep(ti).spikes'];
    end
    [bc,~] = histcounts(spi,bine);
    bc = bc/binw;
    spike_hist(ci,:)=bc/nt; %av across trials
end

%     figure, imagesc(spike_hist)
%     hold on;
%     bin0=interp1(binc,1:nbin,0);
%     plot([bin0,bin0],[0.5, nch+0.5],'--w');
%     xticks(1:2:nbin);
%     xticklabels(binc(1:2:end));
end