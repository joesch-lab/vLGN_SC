function opto_no_other = get_opto_on_edge_no_other(sync_signal, opto_cont, no_before_s, no_after_s, sampling_rate)
%get the time of ON edge of opto flash with no other signal before
[opto_on, ~] = get_signal_on_off_edge(opto_cont);
N=size(sync_signal,2);

%combine opto and vis signal
sigdata = sync_signal | opto_cont;
%find all vis_on edges without opto signal
opto_no_other=[];
onlen=numel(opto_on);
no_signal_before=uint32(no_before_s*sampling_rate);
no_signal_after=uint32(no_after_s*sampling_rate);
for i=1:onlen
    istart=opto_on(i)-no_signal_before;
    if istart<1, continue; end
    iend=opto_on(i)+no_signal_after;   
    if iend>N, continue; end
    if ~any(sigdata(istart:opto_on(i)-1)) && ~any(sync_signal(opto_on(i):iend))
       opto_no_other=[opto_no_other, opto_on(i)];
    end
end

