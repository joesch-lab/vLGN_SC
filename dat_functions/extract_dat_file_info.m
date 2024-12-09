function experiment_info = extract_dat_file_info(datfile)    
    file_info=dir(datfile);
    %get all raw data from dat file     
    ballcam=cell(4,1);
    [camtrig, red_signal, ballcam{1}, ballcam{2}, ballcam{3}, ballcam{4}, optotrig, sample_rate] = read_all_from_Intan_RHD2000_file(datfile);
    experiment_info.sampling_rate = sample_rate;
    
    %get opto values
    if any(optotrig)
        [optoHz, opto_pulse_on_s, opto_pulse_cycle_s, opto_cont] = get_continuous_opto(optotrig, sample_rate);        
    else
        optoHz=NaN; opto_pulse_on_s=0;
        opto_pulse_cycle_s=0; opto_cont=[];
    end
    experiment_info.opto_cont=opto_cont;
    experiment_info.opto_param.optoHz=optoHz;
    experiment_info.opto_param.optoON=opto_pulse_on_s;
    experiment_info.opto_param.opto_cycle=opto_pulse_cycle_s;

    %find cam tigger onset time    
    if datetime(file_info.date)<datetime('Dec-2023')
    experiment_info.cam_trigger =  find(diff([0 camtrig])==1);
    else
        cam_trigger_all =  find(diff([0 camtrig])==1);
        %remove cam triggers before the first and last pause in the camtrig
        %channel
        campause=find(diff(cam_trigger_all)>sample_rate);
        campause1=campause(1)+1;
        campauseN=campause(end);
        camrec=cam_trigger_all >=cam_trigger_all(campause1) & cam_trigger_all <=cam_trigger_all(campauseN);
        experiment_info.cam_trigger =cam_trigger_all(camrec);
%         figure, plot(camtrig)
%         hold on; scatter(cam_trigger_all(camrec), ones(1,sum(camrec)),'*r');
    end

    %red signal  - clean later depending on the stimulus
    experiment_info.red_signal = red_signal; 

    %get loco per sec
    [forward_cm_s, sideways_cm_s, rotation_deg_s, ball_time] = ballcam2speed(ballcam, sample_rate);
    experiment_info.loco.ball_time_s = ball_time;
    experiment_info.loco.forward_cm_s = forward_cm_s;
    experiment_info.loco.sideways_cm_s = sideways_cm_s;
    experiment_info.loco.rotation_deg_s = rotation_deg_s;
end