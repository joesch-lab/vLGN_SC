function opto_yn = is_opto(param_vals_i,pthresh)
opto_yn = param_vals_i.opon_pval<=pthresh || param_vals_i.opoff_pval<=pthresh || ... % opto on or off
    (param_vals_i.flon_pval<=pthresh && param_vals_i.flopon_s_pval<=pthresh) ||...% visual ON responses have been modulated
    (param_vals_i.floff_pval<=pthresh && param_vals_i.flopoff_s_pval<=pthresh) ||...% visual OFF responses have been modulated
    (param_vals_i.flon_z_pval<=pthresh && param_vals_i.flopon_z_pval>pthresh) || ... % visual ON peak became non-significant
    (param_vals_i.floff_z_pval<=pthresh && param_vals_i.flopoff_z_pval>pthresh) ; % visual OFF peak became non-significant
end

