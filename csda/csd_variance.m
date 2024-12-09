function var_cha = csd_variance(CSD,flash_onset_samples)
    sr=20000;        
    nt=size(CSD,2);    
    tafter=0.2*sr;
    tend=min(flash_onset_samples+tafter,nt);
    
%     var_cha=var(CSD(:,flash_onset_samples:tend),[], 2);         
    var_cha=var(CSD,[], 2);         
end