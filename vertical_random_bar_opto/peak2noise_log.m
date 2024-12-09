function logp2n = peak2noise_log(RF, frames_before_spike) 
    %location of the rf as the point with max variance in time
    var_t=var(RF(1:frames_before_spike,:));
    [~,rf_loc]=max(var_t);
    %when the peak happened    
    [~,maxvarT]=max(abs(RF(1:frames_before_spike,rf_loc)));
    peakval=RF(maxvarT,rf_loc);
    var_noise = var(RF(1:10,:),[],"all");    
    peak2noise=peakval^2/var_noise;
    logp2n = 10*log10(peak2noise);
end