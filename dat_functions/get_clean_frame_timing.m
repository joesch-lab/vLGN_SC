function [frame_timing, red_signal_local] = get_clean_frame_timing(red_signal_local, filename)
% removes red frames registered multiple times and inserts missed frames
 
    rftiming=find(red_signal_local);

    %compute the time between neighboring red frames
    drf=rftiming(2:end)-rftiming(1:end-1);
    %intra red frame interval
    ifired=median(drf);
    
    %find red frames with 120ms in between (instead of 83)
    big_gap_duration = 1.5*ifired;
    missed_rf=find(drf>big_gap_duration);     
    %interval in sec*sampling_rate between  frames
    ifi=ifired/5;
    %check if there is a very large gap: several consecutive frames have been missed, results might not be reliable 
    large_missed_rf=find(drf>ifired*2.2); %if more consecutive frames missing issue warning
    if ~isempty(large_missed_rf)
        warning([filename, ': more than one consecutive frame is missing. The results might be wrong.']);    
    end
    if isempty(missed_rf)
        nummissed=0;
    else
        nummissed=length(missed_rf);
    end
    disp([num2str(nummissed),' frames were missed']);
    %insert the missed redfraemes
    while ~isempty(missed_rf)
        % the time where the red frame has to be
        mrf=rftiming(missed_rf)+ifired;
        red_signal_local(mrf)=1;
        rftiming=find(red_signal_local);
        drf=rftiming(2:end)-rftiming(1:end-1);
        missed_rf=find(drf>big_gap_duration); 
    end

    %number of red frames (we do not consider the last red frame and frames after it)
    RFnum=length(rftiming)-1;
    %number of all frames
    allframesnum=RFnum*5;
    %timing of all frames, computed using linear interpolation of time of red
    %frames
    frame_timing=interp1(1:5:allframesnum,rftiming(1:end-1),1:allframesnum,'linear','extrap');
end