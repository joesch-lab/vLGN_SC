function  [optoHz, opto_pulse_on_s, opto_pulse_cycle_s, opto_cont] = get_continuous_opto(optotrig, sampling_rate)
% opto stimulation is done in short pulses (optotrig). This function
% finds contiguous regions when opto stimulation was on
% e.g.: optotrig: 00010101010100001010100000
% will be converted into opto_cont: 00011111111100001111100000
% the function also returns parameters of opto stimulation: Hz, duration of
% on part of the cycle in seconds (opto_pulse_on_s) and duration of the
% whole cycle in secs
%
% O.Symonova, 2022

    opto_on=find(diff([0 optotrig])==1);
    if isempty(opto_on)
        optoHz=NaN;
        opto_pulse_on_s=NaN;
        opto_pulse_cycle_s=NaN;
        opto_cont=NaN;
        return;
    end

    opto_off=find(diff([optotrig 0])==-1);
    opto_pulse_on=median(opto_off-opto_on);
    opto_pulse_on_s=opto_pulse_on/sampling_rate; %how long opto stays on
    opto_pulse_off=median(opto_on(2:end)-opto_off(1:end-1));
    opto_pulse_cycle=opto_pulse_on+opto_pulse_off;
    opto_pulse_cycle_s=opto_pulse_cycle/sampling_rate; %opto cycle : on+off
    optoHz=sampling_rate/opto_pulse_cycle;

    st_id = find(diff(opto_on)>2*opto_pulse_cycle)+1;
    stopto_bout=[opto_on(1),opto_on(st_id)];
    en_id = find(diff(opto_off)>2*opto_pulse_cycle);
    enopto_bout=[opto_off(en_id),opto_off(end)];
    
    opto_cont = zeros(size(optotrig));
    nbouts=length(enopto_bout);
    for i=1:nbouts
        opto_cont(stopto_bout(i):enopto_bout(i))=1;
    end

%     figure, plot(optotrig); hold on;
%     plot(opto_cont); hold on;
%     for i=1:nbouts
%         plot([stopto_bout(i) stopto_bout(i)],[0,1],'--m'); hold on;
%         plot([enopto_bout(i) enopto_bout(i)],[0,1],'--k'); hold on;
%     end

end