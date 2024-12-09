function [flash_on_reps,flash_off_reps, opto_on_reps, opto_off_reps, ...
    flop_on_reps, flop_off_reps] = get_rep_id_offsets(of_data, r_win_s)

%from all offsets of visual&opto flashes find the conditions: flash ON/OFF,
%opto ON/OFF, flash+opto ON/OFF based on the overlap of flash and opto
%intervals

sr=of_data.sr;

opto_dt=r_win_s*sr; %window offset between opto and flash stimulus
r_win=r_win_s*sr;

flash_on_reps =[];
flash_off_reps=[];
opto_on_reps=[];        
opto_off_reps=[];
flop_on_reps=[];        
flop_off_reps=[];

%duration of the bout
tbout=of_data.boutlen_samples;
%start, duration, end of the flash
tf_st=of_data.flash_st;
tf_dur = of_data.flash_dur;
tf_en=tf_st+tf_dur-1;

%start, duration, end of opto
to_st=of_data.flash_st+(of_data.offsets*sr);
to_dur = of_data.opto_dur;
to_en=to_st+to_dur-1;

nofi=length(of_data.offsets);
for ofi=1:nofi 
    % flash on
    sti=tf_st;
    eni=sti+r_win-1;
    if eni<to_st(ofi) %flash response window finished before the start of opto %if ofi>=11
        flash_on_reps=[flash_on_reps,ofi];        
    end

    % flash off
    sti=tf_en+1;
    eni=sti+r_win-1;
    %flash off response window finished before the start of opto or starts 0.2 after opto ends %if ofi<=3
    if eni<=to_st(ofi) || sti> to_en(ofi)+opto_dt
        flash_off_reps=[flash_off_reps,ofi];        
    end

    %opto on
    sti=to_st(ofi);
    eni=sti+r_win-1;
    if eni< tf_st % if ofi<=3   
        opto_on_reps=[opto_on_reps,ofi];         
    end

    %opto off
    sti=to_en(ofi);
    eni=sti+r_win-1;
    if (sti>tf_en+opto_dt || eni<tf_st) && eni<=tbout %ofi>=10                
        opto_off_reps=[opto_off_reps,ofi];        
    end

    %flash on+opto
    sti=tf_st;
    eni=sti+r_win-1;
    if sti>=to_st(ofi) && eni<=to_en(ofi)
        flop_on_reps=[flop_on_reps,ofi];           
    end
    
    %flash off+opto
    sti=tf_en+1;
    eni=sti+r_win-1;    
    if sti>=to_st(ofi) && eni<=to_en(ofi)
        flop_off_reps=[flop_off_reps, ofi];        
    end
end