function vis_yn = is_visual(param_vals_i,pthresh)
vis_yn = param_vals_i.flon_pval<=pthresh || param_vals_i.floff_pval<=pthresh; %flash ON or OFF
