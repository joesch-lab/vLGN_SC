function [all_responses, max_fr_params, spike_count, first_spike, bine, sr] = get_firing_props_per_offset(of_data_file, pval_thresh,arrange_by_intensity)

%for each visual neuron collect flash responses for the intensity: 0->1 and
%1->0 return the array together with the neurong identity and the depth

load(of_data_file);
[dir_of_data,~,~]=fileparts(of_data_file);
[~,exp_name,~]=fileparts(dir_of_data);
params_file=fullfile(dir_of_data,'params_flash_opto_responses_all.mat');
% params_file=fullfile(fullfile(resfolder,exp_name),'params_flash_opto_responses_all.mat');
load(params_file);

spike_count=[];
first_spike=[];
all_responses=[];
%params: [ vis_max_fr_off, vis_tmax_fr_off, vis_width_peak_off, vis_z_pval_off,...
%  opto_max_fr_off, opto_tmax_fr_off, opto_width_peak_off, opto_z_pval_off,...
%      opto_long_max_fr_off, opto_long_tmax_fr_off, opto_long_width_peak_off, opto_long_z_pval_off];
max_fr_params=[];


%params for data extraction
sr=of_data.sr;
r_win_s=0.2; %window to compute flash on, flash off, opto on, opto off responses
tbefore_s=r_win_s/2;

tf_st=of_data.flash_st;
tf_dur = of_data.flash_dur;
tf_en=tf_st+tf_dur-1;
tf_st_s=tf_st/sr;
tf_en_s=tf_en/sr;

%start, duration, end of opto
to_st=of_data.flash_st+(of_data.offsets*sr);
to_dur = of_data.opto_dur;
to_en=to_st+to_dur-1;
to_st_s=to_st/sr;
to_en_s=to_en/sr;


bin_dur_s=of_data.histograms.binlen; %bin duration to compute responses
bine=-tbefore_s:bin_dur_s:r_win_s;
nbinbefore = tbefore_s/bin_dur_s;
nbintafter = r_win_s/bin_dur_s;

nofi=size(of_data.offsets,2);
%collect reps where there is no overlap with opto
[flash_on_reps,flash_off_reps, opto_on_reps, opto_off_reps, ...
    flop_on_reps, flop_off_reps] = get_rep_id_offsets(of_data, r_win_s);

ncl=length(of_data.spikes);
for ci=1:ncl
    %select only visual and opto neurons
    if ~(is_visual(params_vals(ci),pval_thresh) && is_opto(params_vals(ci),pval_thresh))
        continue;
    end

    %collect spikes
    first_spikes_flash_on=[];
    first_spikes_flash_off=[];
    first_spikes_flop_on=[];
    first_spikes_flop_off=[];
    first_spikes_opto_on=[];
    first_spikes_opto_off=[];

    %collect histograms from ofdata
    flon_hist_i=[];
    floff_hist_i=[];
    flopon_hist_i=[];
    flopoff_hist_i=[];
    opon_hist_i=[];
    opoff_hist_i=[];

    for ofi=1:nofi
        nrep_i=size(of_data.spikes(ci).offset(ofi).rep,2);
        if ismember(ofi,flash_on_reps)
            sti=of_data.histograms.bin_flash_st-nbinbefore+1;
            eni = of_data.histograms.bin_flash_st+nbintafter;
            flon_hist_i=[flon_hist_i; reshape(of_data.histograms.hist(ci,ofi,sti:eni),1,[])];
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp;
                spi=spi/sr - tf_st_s;
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                first_spikes_flash_on=[first_spikes_flash_on, min(spi(spi>0))];
            end
        end
        if ismember(ofi,flash_off_reps)
            sti=of_data.histograms.bin_flash_en-nbinbefore+1;
            eni = of_data.histograms.bin_flash_en+nbintafter;
            floff_hist_i=[floff_hist_i; reshape(of_data.histograms.hist(ci,ofi,sti:eni),1,[])];
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp;
                spi=spi/sr - tf_en_s;
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                first_spikes_flash_off=[first_spikes_flash_off, min(spi(spi>0))];
            end
        end
        if ismember(ofi,flop_on_reps)
            sti=of_data.histograms.bin_flash_st-nbinbefore+1;
            eni = of_data.histograms.bin_flash_st+nbintafter;
            flopon_hist_i=[flopon_hist_i; reshape(of_data.histograms.hist(ci,ofi,sti:eni),1,[])];
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp;
                spi=spi/sr - tf_st_s;
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                first_spikes_flop_on=[first_spikes_flop_on, min(spi(spi>0))];
            end
        end
        if ismember(ofi,flop_off_reps)
            sti=of_data.histograms.bin_flash_en-nbinbefore+1;
            eni = of_data.histograms.bin_flash_en+nbintafter;
            flopoff_hist_i=[flopoff_hist_i; reshape(of_data.histograms.hist(ci,ofi,sti:eni),1,[])];
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp;
                spi=spi/sr - tf_en_s;
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                first_spikes_flop_off=[first_spikes_flop_off, min(spi(spi>0))];
            end
        end
        if ismember(ofi,opto_on_reps)
            sti=round(of_data.histograms.bin_opto_st(ofi))-nbinbefore+1;
            eni = round(of_data.histograms.bin_opto_st(ofi))+nbintafter;
            opon_hist_i=[opon_hist_i; reshape(of_data.histograms.hist(ci,ofi,sti:eni),1,[])];
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp;
                spi=spi/sr - to_st_s(ofi);
                spi_bl=spi(spi>=-r_win_s & spi<0);
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                first_spikes_opto_on=[first_spikes_opto_on, min(spi(spi>0))];
            end
        end
        if ismember(ofi,opto_off_reps)
            sti=round(of_data.histograms.bin_opto_en(ofi))-nbinbefore+1;
            eni = round(of_data.histograms.bin_opto_en(ofi))+nbintafter;
            opoff_hist_i=[opoff_hist_i; reshape(of_data.histograms.hist(ci,ofi,sti:eni),1,[])];
            for rii=1:nrep_i
                spi=of_data.spikes(ci).offset(ofi).rep(rii).sp;
                spi=spi/sr - to_en_s(ofi);
                spi=spi((spi>-tbefore_s & spi<=r_win_s));
                first_spikes_opto_off=[first_spikes_opto_off, min(spi(spi>0))];
            end
        end
    end

    %mean histogram across reps
    flash_on_hist = mean(flon_hist_i,1);
    flash_off_hist = mean(floff_hist_i,1);
    opto_on_hist = mean(opon_hist_i,1);
    opto_off_hist = mean(opoff_hist_i,1);
    flop_on_hist = mean(flopon_hist_i,1);
    flop_off_hist = mean(flopoff_hist_i,1);

    bin0=floor(interp1(bine(1:end-1),1:length(bine)-1,0));
    bin_01=floor(interp1(bine(1:end-1),1:length(bine)-1,0.1));

    spike_count_i=[sum(flash_on_hist(bin0:end))*bin_dur_s,...
        sum(flash_off_hist(bin0:end))*bin_dur_s,...
        sum(flop_on_hist(bin0:end))*bin_dur_s,...
        sum(flop_off_hist(bin0:end))*bin_dur_s,...
        sum(opto_on_hist(bin0:end))*bin_dur_s,...
        sum(opto_off_hist(bin0:end))*bin_dur_s,...
        sum(opto_on_hist(1:bin0-1))*bin_dur_s,... %0.1s of baseline before opto
        sum(opto_on_hist(bin0:bin_01-1))*bin_dur_s]; %0.1s of of opto

    %arrange in this shape: [white_flash, black_flash, wf+opto, bf+opto, opto on, opto off]
    if arrange_by_intensity && of_data.flashI<10 %black flash on-flash evokes off-responses
        %black flash, get on response
        all_responses=[all_responses;[ci,of_data.cluster_depth(ci),flash_off_hist,flash_on_hist, flop_off_hist, flop_on_hist, opto_on_hist, opto_off_hist]];
        parami=[params_vals(ci).flon_z_frpeak_maxfr, params_vals(ci).flon_z_frpeak_t, params_vals(ci).flon_z_frpeak_w, params_vals(ci).flon_z_pval,...
            params_vals(ci).flopon_z_frpeak_maxfr, params_vals(ci).flopon_z_frpeak_t, params_vals(ci).flopon_z_frpeak_w, params_vals(ci).flopon_z_pval,...
            params_vals(ci).flopon_long_z_frpeak_maxfr, params_vals(ci).flopon_long_z_frpeak_t, params_vals(ci).flopon_long_z_frpeak_w, params_vals(ci).flopon_long_z_pval];
        spike_count=[spike_count;[ci,of_data.cluster_depth(ci),spike_count_i(2),spike_count_i(1), spike_count_i(4), spike_count_i(3),spike_count_i(5:end)]];
        first_spike = [first_spike; [ci,of_data.cluster_depth(ci),median(first_spikes_flash_off),median(first_spikes_flash_on), median(first_spikes_flop_off), median(first_spikes_flop_on), median(first_spikes_opto_on),median(first_spikes_opto_off)]];
    else %white flash: on-flash for on-responses
        all_responses=[all_responses;[ci,of_data.cluster_depth(ci),flash_on_hist, flash_off_hist,flop_on_hist, flop_off_hist,opto_on_hist, opto_off_hist]];
        parami=[params_vals(ci).floff_z_frpeak_maxfr, params_vals(ci).floff_z_frpeak_t, params_vals(ci).floff_z_frpeak_w, params_vals(ci).floff_z_pval,...
            params_vals(ci).flopoff_z_frpeak_maxfr, params_vals(ci).flopoff_z_frpeak_t, params_vals(ci).flopoff_z_frpeak_w, params_vals(ci).flopoff_z_pval,...
            params_vals(ci).flopoff_long_z_frpeak_maxfr, params_vals(ci).flopoff_long_z_frpeak_t, params_vals(ci).flopoff_long_z_frpeak_w, params_vals(ci).flopoff_long_z_pval];
        spike_count=[spike_count;[ci,of_data.cluster_depth(ci),spike_count_i]];
        first_spike = [first_spike; [ci,of_data.cluster_depth(ci),median(first_spikes_flash_on), median(first_spikes_flash_off),median(first_spikes_flop_on), median(first_spikes_flop_off), median(first_spikes_opto_on),median(first_spikes_opto_off)]];
    end    
    max_fr_params=cat(1,max_fr_params,parami);
end
end

