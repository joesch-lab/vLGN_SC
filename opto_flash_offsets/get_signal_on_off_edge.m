function [on_e, off_e] = get_signal_on_off_edge(signal, windt)
    %redframe with the next 3 frames
    if ~exist('windt','var')
        windt=uint16(0.0167*20000);
    end
    
    %moving sum, include windth values on the rigt   
    signal_sum=movsum(signal,[0,windt]);
    flashvals=signal_sum>0;
    %difference between neighboring values
    signal_diff=diff([flashvals, 0]);
    off_e=find(signal_diff==-1);
    
    %moving sum, include windt values on the left
    signal_sum=movsum(signal,[windt,0]);
    flashvals=signal_sum>0;
    %difference between neighboring values
    signal_diff=diff([0, flashvals]);
    on_e=find(signal_diff==1);
    % on_signal=zeros(size(signal));
    % off_signal=zeros(size(signal));
    % on_signal(on_e)=1;
    % off_signal(off_e)=1;
    % figure, plot(signal); hold on; plot(0.5*on_signal); hold on; plot(0.5*off_signal);
end