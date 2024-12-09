function [spike_hist, spike_count, spike_count_before, binc] = get_spikes_on_times_perrep(tpts,spike_data, sr, tbefore, tafter, binw0, spike_count_win)
% compute spike histogram and counts around specific time points
% params:
% tpts: time points around which to compute histograms
% spike_data: data structure with spikes and unit information
% sr: sampling rate of spikes (20K for silicone probe, 50K for neuropixels)
% tbefore: time in secs to consider before event onset
% tafter: time in secs  to consider after event onset, for histogram computation
% binw0: width of the bin
% spike_count_win: interval in secs for spike count before and after event onset
% output:
% spike_hist: PSTH for each cluster, reps not averaged
% spike_count: spike count per each cluster [0..spike_count_win],
% spike_count_before: spike count per each cluster [-spike_count_win..0],
% binc: center of the bins in spike histograms
%
% O.Symonova, 2024

nch=length(spike_data);
nrep=length(tpts);
spikes_raw={};
for ci=1:nch
    if abs(sr-20000)<1
        spikest=double(spike_data(ci).spike_timing)/sr;
    else
        spikest=double(spike_data(ci).spike_timing);
    end
    for ti=1:nrep
        tst=tpts(ti)-tbefore;
        ten=tpts(ti)+tafter;
        sp_sel=(spikest(spikest>=tst & spikest<=ten))-tpts(ti);
        spikes_raw(ci).rep(ti).spikes=sp_sel;
    end
end

%     for ci=1:nch
%         nsp=0;
%         for ti=1:nrep
%             nsp=nsp+numel(spikes_raw(ci).rep(ti).spikes);
%         end
%         if nsp==0
%             continue;
%         end
%         f=figure;
%         for ti=1:nrep
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

spike_hist=zeros(nrep, nch,nbin);
spike_count=zeros(nrep, nch);
spike_count_before=zeros(nrep, nch);
for ti=1:nrep
    for ci=1:nch
        spi=spikes_raw(ci).rep(ti).spikes';
        [bc,~] = histcounts(spi,bine);
        spike_hist(ti,ci,:)= bc/binw;
        spike_count(ti,ci)=sum(spi>=0 & spi<=spike_count_win);
        spike_count_before(ti,ci)=sum(spi>=-spike_count_win & spi<=0);
    end
end

end