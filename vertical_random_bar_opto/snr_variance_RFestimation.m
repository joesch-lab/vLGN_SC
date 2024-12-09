function snr_variance = snr_variance_RFestimation(RF, frames_before_spike, twin_around_max) 
        var_t=var(RF(1:frames_before_spike,:),[],2);
        [~,t_maxvar]=max(var_t);
        tmin=max(1,t_maxvar-twin_around_max);
        tmax=max(1,t_maxvar+twin_around_max);
        %variance at the max t
        varsignal=var(RF(tmin:tmax,:),[],"all");
        %variance of noise
        varnoise=var(RF([1:tmin-1,tmax+1:frames_before_spike],:),[],"all");
        snr_variance =  varsignal/varnoise;
end