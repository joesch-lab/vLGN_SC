function vis_no_other = get_vis_on_edge_no_other(sync_signal, opto_cont, no_before_s, no_after_s, sampling_rate)
%get the time of ON edge of opto flash with no other signal before
[vis_on, ~] = get_signal_on_off_edge(sync_signal);
N=size(sync_signal,2);

%combine opto and vis signal
sigdata = sync_signal | opto_cont;
%find all vis_on edges without opto signal
vis_no_other=[];
onlen=numel(vis_on);
no_signal_before=uint32(no_before_s*sampling_rate);
no_signal_after=uint32(no_after_s*sampling_rate);
for i=1:onlen
    istart=vis_on(i)-no_signal_before;
    if istart<1, continue; end;
    iend=vis_on(i)+no_signal_after;   
    if iend>N, continue; end;
    if ~any(sigdata(istart:vis_on(i)-1)) && ~any(opto_cont(vis_on(i):iend))
        vis_no_other=[vis_no_other, vis_on(i)];
    end
end

