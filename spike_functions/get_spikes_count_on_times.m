function [spike_count, cluster_depth] = get_spikes_count_on_times(tpts,spike_data, t_secs)
% get spike count per unit in a window [0,t_secs] after event onset
% params:
% tpts: time points around which to compute spike count
% spike_data: data structure with spikes and unit information
% spike_count_win: interval in secs for spike count after event onset
%
% O.Symonova, 2024

    nch=length(spike_data);
    nt=length(tpts);
    spike_count=zeros(nch,nt);
    cluster_depth=zeros(nch,1);
    if isa(spike_data(1).spike_timing(1), 'uint64')
        istntdata=false;
    else
        istntdata=true;
    end

    for ci=1:nch         
        for ti=1:nt
            if ~istntdata
                spikest=double(spike_data(ci).spike_timing)/20000;
            else
                spikest=double(spike_data(ci).spike_timing);
            end
            tst=tpts(ti);
            ten=tpts(ti)+t_secs;            
            spike_count(ci,ti)=sum(spikest>=tst & spikest<=ten);
            cluster_depth(ci)=spike_data(ci).SC_depth;
        end
    end