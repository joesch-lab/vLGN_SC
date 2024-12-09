function var_t = CSDA_response_timing(csd_arr, onset, before_s, after_s, sr)
before_sr=before_s*sr;
after_sr=after_s*sr;
st=onset-before_sr;
en = min(onset + after_sr,size(csd_arr,2));
var_t  = var(csd_arr,'omitnan');
var_before_onset=mean(var_t(st:onset));
var_t=var_t-var_before_onset;
var_t=var_t(st:en);
end

