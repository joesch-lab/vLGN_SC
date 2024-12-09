function compute_response_params(of_data_file, resfolder)
% for each cluster computer whether if is responsive to the visual & opto
% flashes

load(of_data_file);
params_file=fullfile(resfolder,'params_flash_opto_responses_all.mat');

zeta_params={};
zeta_params.intResampNum = 1000; %~50 random resamplings should give us a good enough idea if this cell is responsive, but if it's close to 0.05, we should increase this #. Generally speaking, more is better, so let's put 100 here.
zeta_params.intPlot = 0;%what do we want to plot?(0=nothing, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot)
zeta_params.vecRestrictRange = [0 inf];%do we want to restrict the peak detection to for example the time during stimulus? Then put [0 1] here.
zeta_params.boolDirectQuantile = false;%if true; uses the empirical null distribution rather than the Gumbel approximation. Note that in this case the accuracy of your p-value is limited by the # of resamplings
zeta_params.dblJitterSize = 1; %scalar value, sets the temporal jitter window relative to dblUseMaxDur (default: 2). Note that a value of "2" therefore means that the last data point that can be used is at t=last_onset + dblUseMaxDur*3
zeta_params.boolStitch = false;

%params for data extraction
sr=of_data.sr;
opto_dt=0.2*sr; %window offset between opto and flash stimulus
r_win_s=0.2; %window to compute flash on, flash off, opto on, opto off responses
tb=0.025; %small window to exclude computation before/after flash/opto
r_win=r_win_s*sr;

%duration of the bout
tbout=of_data.boutlen_samples;
%start, duration, end of the flash
tf_st=of_data.flash_st;
tf_dur = of_data.flash_dur;
tf_en=tf_st+tf_dur-1;

%start, duration, end of opto
to_st=of_data.flash_st+(of_data.offsets*sr);
to_dur = of_data.opto_dur;
to_en=min(tbout,to_st+to_dur-1);

%number of different offsets
nofi=size(of_data.offsets,2);

%count total number of reps
nreps_all=0;
for ofi=1:13
    nreps_all=nreps_all+size(of_data.spikes(1).offset(ofi).rep,2);
end

%arrays to save full flash and opto duration
flash_on_full=zeros(nreps_all,tbout);
opto_on_full=zeros(nreps_all,tbout);
% baseline, flash-on, flash-off, opto-on, opto-off, flash on + opto, ...
% flash off + opto, flash on + opto (long effect), flash off + opto (long effect)
b_array=zeros(nreps_all,tbout);
flon_array=zeros(nreps_all,tbout);
floff_array=zeros(nreps_all,tbout);
opon_array=zeros(nreps_all,tbout);
opoff_array=zeros(nreps_all,tbout);
flopon_array=zeros(nreps_all,tbout);
flopoff_array=zeros(nreps_all,tbout);
flopon_long_array=zeros(nreps_all,tbout);
flopoff_long_array=zeros(nreps_all,tbout);

% baseline: all lines before start of the stimulus
% flash on: [flash_start : flash_start + r_win] (in the last 3 offsets before opto, max distanse from opto)
% flash off: [flash_end : flash_end + r_win] (in the first 3 offsets after opto, max distanse from opto)
% opto-on: [opto_start : + r_win] (before flash, first 3 offsets)
% opto-off: [opto_off : opto_off + r_win] (last 4 offsets, after flash + exlude_time)
% flash on + opto: [flash_start: flash_start + r_win] (during 5,6,7, offsets, when flash on and opto overlap)
% flash off + opto: [flash_end: flash_end + r_win] (during 9,10,11 offsets, when flash off and opto overlap)
% flopon_long: [flash_start: flash_start + r_win] (flash-on resp is in the window [0, 0.2]s after opto release)
% flopoff_long: [flash_start: flash_end + r_win] (flash-off resp is in the window [0, 0.2]s after opto release)

ri=1;
for ofi=1:nofi  
    nrep_i=size(of_data.spikes(1).offset(ofi).rep,2);

    %find the start off the stimuli during current offset
    stim_start=max(1,min(tf_st,to_st(ofi))-tb*sr);
    %baseline is before any stimulus        
    b_array(ri:ri+nrep_i-1,1:stim_start)=1;    
    
    % flash on
    sti=tf_st;
    eni=sti+r_win-1;   
    if eni<to_st(ofi) %flash response window finished before the start of opto %if ofi>=11
        flon_array(ri:ri+nrep_i-1,sti:eni)=1;        
    end

    %full flash on: for visualization
    flash_on_full(:,tf_st:tf_en)=1;

    % flash off
    sti=tf_en+1;
    eni=sti+r_win-1;
    %flash off response window finished before the start of opto or starts 0.2 after opto ends %if ofi<=3
    if eni<=to_st(ofi) || sti> to_en(ofi)+opto_dt 
        floff_array(ri:ri+nrep_i-1,sti:eni)=1;
    end

    %opto on
    sti=to_st(ofi);
    eni=sti+r_win-1;
    if eni< tf_st % if ofi<=3                
        opon_array(ri:ri+nrep_i-1,sti:eni)=1;
    end

    %full opto on: for visualization        
    sti=to_st(ofi);
    eni=to_en(ofi);       
    opto_on_full(ri:ri+nrep_i-1,sti:eni)=1;    

    %opto off
    sti=to_en(ofi);
    eni=min(tbout,sti+r_win-1);
    if sti>tf_en+ opto_dt || eni<tf_st %ofi>=10                
        opoff_array(ri:ri+nrep_i-1,sti:eni)=1;
    end

    %flash on+opto
    sti=tf_st;
    eni=sti+r_win-1;        
    if sti>=to_st(ofi) && eni<=to_en(ofi)        
        flopon_array(ri:ri+nrep_i-1,sti:eni)=1;    
    end

    %flash off+opto
    sti=tf_en+1;
    eni=sti+r_win-1;
    if sti>=to_st(ofi) && eni<=to_en(ofi)
        flopoff_array(ri:ri+nrep_i-1,sti:eni)=1;
    end


%     %flash on+opto
%     if ofi>4 && ofi<8
%         sti=tf_st;
%         eni=sti+r_win-1;        
%         flopon_array(ri:ri+nrep_i-1,sti:eni)=1;    
%     end
% 
%     %flash off+opto
%     if ofi>8 && ofi<12
%         sti=tf_en;
%         eni=sti+r_win-1;
%         flopoff_array(ri:ri+nrep_i-1,sti:eni)=1;
%     end

    %flash on+ opto (long effect) 
    %flash response starts after opto and ends within opto_dt after opto release
    sti=tf_st;
    eni=sti+r_win-1;       
    if sti>=to_en(ofi) && eni<=to_en(ofi)+opto_dt  
        flopon_long_array(ri:ri+nrep_i-1,sti:eni)=1;        
    end

    %flash off+opto (long effect)
    %flash off response starts after opto and ends within opto_dt after opto release
    sti=tf_en+1;
    eni=sti+r_win-1;
    if sti>=to_en(ofi) && eni<=to_en(ofi)+opto_dt  
        flopoff_long_array(ri:ri+nrep_i-1,sti:eni)=1;        
    end
    ri=ri+nrep_i;
end

%% visualize time intervals for computing flash/opto responses
imvis = zeros(size(b_array));
imvis(b_array>0)=1;
imvis(flash_on_full>0)=2;
imvis(opto_on_full>0)=9;
imvis(flon_array>0)=3;
imvis(flopon_array>0)=4;
imvis(flopon_long_array>0)=5;

imvis(floff_array>0)=6;
imvis(flopoff_array>0)=7;
imvis(flopoff_long_array>0)=8;

imvis(opon_array>0)=10;
imvis(opoff_array>0)=11;



figure, imagesc(imvis);
set(gca,'YDir','normal');
cmap = [255,255,255;... %all
    204,204,204;...%baseline
    204,255,204;... % full flash
    111,255,111; ... %flash on
    119,216,119;... %flash on +opto
    98,178,98;... %flash on + opto long   
    133,234,255;... % flash off
    0,195,255;... %flash off + opto 
    0,148,178;... %flash off + opto long    
    204,204,255;... %full opto
    209,125,249;... %opto on
    144,86,172;... %opto off    
    ]/255;

colormap(cmap);
hcb = colorbar;
hcb.Ticks=linspace(0.5,10.5,12);
hcb.TickLength = 0; 
hcb.TickLabels=[{'all'},{'baseline'},{'flash'},{'flash on'},{'flash on + opto'},...
{'flash on + opto s.'}, {'flash off'},{'flash off + opto'},{'flash off + opto s.'},...
{'opto'},{'opto on'},{'opto off'}];


%% get the onset time of each rep
%flash on
rinds=find(any(flon_array,2));
flash_on_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    flash_on_st(ri)=((rinds(ri)-1)*tbout+find(flon_array(rinds(ri),:),1,"first"));
end
%flash off
rinds=find(any(floff_array,2));
flash_off_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    flash_off_st(ri)=((rinds(ri)-1)*tbout+find(floff_array(rinds(ri),:),1,"first"));
end
%opto on
rinds=find(any(opon_array,2));
opto_on_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    opto_on_st(ri)=((rinds(ri)-1)*tbout+find(opon_array(rinds(ri),:),1,"first"));
end
%opto off
rinds=find(any(opoff_array,2));
opto_off_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    opto_off_st(ri)=((rinds(ri)-1)*tbout+find(opoff_array(rinds(ri),:),1,"first"));
end
%flash on + opto
rinds=find(any(flopon_array,2));
flop_on_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    flop_on_st(ri)=((rinds(ri)-1)*tbout+find(flopon_array(rinds(ri),:),1,"first"));
end
%flash on + long opto
rinds=find(any(flopon_long_array,2));
flop_on_long_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    flop_on_long_st(ri)=((rinds(ri)-1)*tbout+find(flopon_long_array(rinds(ri),:),1,"first"));
end
%flash off + opto
rinds=find(any(flopoff_array,2));
flop_off_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    flop_off_st(ri)=((rinds(ri)-1)*tbout+find(flopoff_array(rinds(ri),:),1,"first"));
end
%flash off + opto
rinds=find(any(flopoff_long_array,2));
flop_off_long_st=zeros(length(rinds),1);
for ri=1:length(rinds)
    flop_off_long_st(ri)=((rinds(ri)-1)*tbout+find(flopoff_long_array(rinds(ri),:),1,"first"));
end



% for all reps compute mean fr to the following conditions:
% baseline, flash ON, flash OFF, opto ON, opto OFF,... 
% flash ON + opto, flash OFF+ opto, flash ON + opto sustained, flash OFF+ opto sustained;  
% compare the difference of means with permutation tests 
% additionally compute zeta-test to capture neurons with similar firing
% rate but event-locked firing events


%duration or reps to convert spike counts to firing rates
flon_dur_reps=sum(flon_array,2)/sr;
floff_dur_reps=sum(floff_array,2)/sr;
opon_dur_reps=sum(opon_array,2)/sr;
opoff_dur_reps=sum(opoff_array,2)/sr;
flopon_dur_reps=sum(flopon_array,2)/sr;
flopoff_dur_reps=sum(flopoff_array,2)/sr;
flopon_long_dur_reps=sum(flopon_long_array,2)/sr;
flopoff_long_dur_reps=sum(flopoff_long_array,2)/sr;

%get samples from the baseline to be used in the randomization tests
n_tests=1000;
[row_bl_off, sti_bl_off]=get_bl_samples(b_array,n_tests,r_win);

%array to store all spikes
spikes=zeros(nreps_all,tbout);

params_vals={};
ncl=length(of_data.spikes);
for ci=1:ncl    
%     disp(ci);
    %collect spikes
    spikes(:)=0;
    ri=1;
    for ofi=1:nofi        
        nrep_i=size(of_data.spikes(ci).offset(ofi).rep,2);        
        for rii=1:nrep_i
            spikes(ri,of_data.spikes(ci).offset(ofi).rep(rii).sp)=1;
            ri=ri+1;
        end
    end
    spikes=spikes(:,1:tbout);
   
    %% baseline fr, randomly sampled and average    
    bl_fr_samples = get_fr_samples(spikes, r_win, row_bl_off, sti_bl_off, sr);
    params_vals(ci).bl_meanfr=mean(bl_fr_samples,"omitnan");  

    %% flash on
    %zeta test flash on
    [rowi,coli]=find(spikes.*flon_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
    
%     [dblZetaP,vecLatencies,sZETA,sRate] = getZeta(ti,flash_on_st/sr,r_win/sr); 
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,flash_on_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch);

    params_vals(ci).flon_z_pval=dblZetaP;
    params_vals(ci).flon_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).flon_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).flon_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).flon_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).flon_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for flash on vs baseline
    %firing rates for reps: flash ON    
    flash_on_vals=sum(spikes.*flon_array,2);
    flash_on_vals=flash_on_vals(flon_dur_reps>0);
    flash_on_vals=flash_on_vals./flon_dur_reps(flon_dur_reps>0); %firing rate per sec
    
    [pval_s, mean_diff, effsize] = permutationTest_with_subsampling(flash_on_vals, bl_fr_samples,n_tests);
    params_vals(ci).flon_meanfr=mean(flash_on_vals,"omitnan");        
    params_vals(ci).flon_s_pval=pval_s;
    params_vals(ci).flon_s_effsize=effsize; 
    params_vals(ci).flon_s_obsdiff=mean_diff;    
    params_vals(ci).flon_pval=min(params_vals(ci).flon_s_pval,params_vals(ci).flon_z_pval);
    
    %% flash off
    %zeta test flash off
    [rowi,coli]=find(spikes.*floff_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);      
    
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,flash_off_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch);
    params_vals(ci).floff_z_pval=dblZetaP;
    params_vals(ci).floff_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).floff_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).floff_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).floff_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).floff_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for flash off vs baseline
    %firing rates for flash OFF    
    flash_off_vals=sum(spikes.*floff_array,2);
    flash_off_vals=flash_off_vals(floff_dur_reps>0);
    flash_off_vals=flash_off_vals./floff_dur_reps(floff_dur_reps>0); %firing rate per sec

    [pval_s, mean_diff, effsize] = permutationTest_with_subsampling(flash_off_vals, bl_fr_samples,n_tests);    
    params_vals(ci).floff_meanfr=mean(flash_off_vals,"omitnan");        
    params_vals(ci).floff_s_pval=pval_s;
    params_vals(ci).floff_s_effsize=effsize; 
    params_vals(ci).floff_s_obsdiff=mean_diff;    
    params_vals(ci).floff_pval=min(params_vals(ci).floff_s_pval,params_vals(ci).floff_z_pval);

    %% opto on
    %zeta test opto on
    [rowi,coli]=find(spikes.*opon_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
        
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,opto_on_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch);
    params_vals(ci).opon_z_pval=dblZetaP;
    params_vals(ci).opon_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).opon_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).opon_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).opon_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).opon_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for opto on vs baseline
    %firing rates for opto ON  
    opto_on_vals=sum(spikes.*opon_array,2);
    opto_on_vals=opto_on_vals(opon_dur_reps>0);
    opto_on_vals=opto_on_vals./opon_dur_reps(opon_dur_reps>0); %firing rate per sec

    [pval_s, mean_diff, effsize]  = permutationTest_with_subsampling(opto_on_vals, bl_fr_samples,n_tests);  
    params_vals(ci).opon_meanfr=mean(opto_on_vals,"omitnan");        
    params_vals(ci).opon_s_pval=pval_s;
    params_vals(ci).opon_s_effsize=effsize; 
    params_vals(ci).opon_s_obsdiff=mean_diff;    
    params_vals(ci).opon_pval=min(params_vals(ci).opon_s_pval,params_vals(ci).opon_z_pval);

    %% opto off
    %zeta test opto off
    [rowi,coli]=find(spikes.*opoff_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
    
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,opto_off_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch); 
    params_vals(ci).opoff_z_pval=dblZetaP;
    params_vals(ci).opoff_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).opoff_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).opoff_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).opoff_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).opoff_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for opto off vs baseline
    %firing rates for opto OFF    
    opto_off_vals=sum(spikes.*opoff_array,2);
    opto_off_vals=opto_off_vals(opoff_dur_reps>0);
    opto_off_vals=opto_off_vals./opoff_dur_reps(opoff_dur_reps>0); %firing rate per sec    

    [pval_s, mean_diff, effsize]  = permutationTest_with_subsampling(opto_off_vals, bl_fr_samples,n_tests);
    params_vals(ci).opoff_meanfr=mean(opto_off_vals,"omitnan");        
    params_vals(ci).opoff_s_pval=pval_s;
    params_vals(ci).opoff_s_effsize=effsize; 
    params_vals(ci).opoff_s_obsdiff=mean_diff;    
    params_vals(ci).opoff_pval=min(params_vals(ci).opoff_s_pval,params_vals(ci).opoff_z_pval);

    %% flash on+opto
    %zeta test flash on+opto
    [rowi,coli]=find(spikes.*flopon_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
    
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,flop_on_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch);
    params_vals(ci).flopon_z_pval=dblZetaP;
    params_vals(ci).flopon_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).flopon_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).flopon_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).flopon_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).flopon_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for flash_on +opto vs flash_on
    %fr for flash ON + opto per each rep
    flash_on_opto_vals=sum(spikes.*flopon_array,2);
    flash_on_opto_vals=flash_on_opto_vals(flopon_dur_reps>0);
    flash_on_opto_vals=flash_on_opto_vals./flopon_dur_reps(flopon_dur_reps>0); %firing rate per sec
        
    %shuffle: flash on vs. flash on+opto
    [pval_s, mean_diff, effsize]  = permutationTest_with_subsampling(flash_on_vals, flash_on_opto_vals,n_tests);
    params_vals(ci).flopon_meanfr=mean(flash_on_opto_vals,"omitnan");
    params_vals(ci).flopon_s_pval=pval_s;    
    params_vals(ci).flopon_s_effsize=effsize; 
    params_vals(ci).flopon_s_obsdiff=mean_diff;

    %% flash on+opto sustained effect
    %zeta test flash on+ long opto
    [rowi,coli]=find(spikes.*flopon_long_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
    
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,flop_on_long_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch); 
    params_vals(ci).flopon_long_z_pval=dblZetaP;
    params_vals(ci).flopon_long_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).flopon_long_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).flopon_long_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).flopon_long_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).flopon_long_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for flash_on +long opto vs flash_on
    %fr for flash ON + long opto per each rep
    flash_on_opto_long_vals=sum(spikes.*flopon_long_array,2);
    flash_on_opto_long_vals=flash_on_opto_long_vals(flopon_long_dur_reps>0);
    flash_on_opto_long_vals=flash_on_opto_long_vals./flopon_long_dur_reps(flopon_long_dur_reps>0); %firing rate per sec
        
    %shuffle: flash on vs. flash on+long opto
    [pval_s, mean_diff, effsize]  = permutationTest_with_subsampling(flash_on_vals, flash_on_opto_long_vals,n_tests);
    params_vals(ci).flopon_long_meanfr=mean(flash_on_opto_long_vals,"omitnan");
    params_vals(ci).flopon_long_s_pval=pval_s;    
    params_vals(ci).flopon_long_s_effsize=effsize; 
    params_vals(ci).flopon_long_s_obsdiff=mean_diff;

    %% flash off+opto
    %zeta test flash off+opto
    [rowi,coli]=find(spikes.*flopoff_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
    
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,flop_off_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch);
    params_vals(ci).flopoff_z_pval=dblZetaP;
    params_vals(ci).flopoff_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).flopoff_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).flopoff_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).flopoff_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).flopoff_z_frpeak_w = sRate.dblPeakWidth;
   
    %fr for flash OFF + opto per each rep    
    flash_off_opto_vals=sum(spikes.*flopoff_array,2);
    flash_off_opto_vals=flash_off_opto_vals(flopoff_dur_reps>0);
    flash_off_opto_vals=flash_off_opto_vals./flopoff_dur_reps(flopoff_dur_reps>0); %firing rate per sec
     
    %shuffle: flash off vs flash off+opto
    [pval_s, mean_diff, effsize]  = permutationTest_with_subsampling(flash_off_vals, flash_off_opto_vals,n_tests);           
    params_vals(ci).flopoff_meanfr=mean(flash_off_opto_vals,"omitnan");
    params_vals(ci).flopoff_s_pval=pval_s;    
    params_vals(ci).flopoff_s_effsize=effsize; 
    params_vals(ci).flopoff_s_obsdiff=mean_diff; 

    %% flash off+opto sustained effect
    %zeta test flash off+ long opto
    [rowi,coli]=find(spikes.*flopoff_long_array);
    ti=sort(((rowi-1)*tbout+coli)/sr);  
    
    [dblZetaP,sZETA,sRate,sLatencies] = zetatest(ti,flop_off_long_st/sr,r_win/sr,...
        zeta_params.intResampNum,zeta_params.intPlot,zeta_params.vecRestrictRange,zeta_params.boolDirectQuantile,zeta_params.dblJitterSize,zeta_params.boolStitch);
    params_vals(ci).flopoff_long_z_pval=dblZetaP;
    params_vals(ci).flopoff_long_z_peak=sZETA.dblZETA; % deviation from uniform
    params_vals(ci).flopoff_long_z_tpeak=sZETA.dblZetaT; %time of the peak
    params_vals(ci).flopoff_long_z_frpeak_t = sRate.dblPeakTime; 
    params_vals(ci).flopoff_long_z_frpeak_maxfr = sRate.dblPeakRate; 
    params_vals(ci).flopoff_long_z_frpeak_w = sRate.dblPeakWidth;

    %shuffle test for flash_on +long opto vs flash_on
    %fr for flash ON + long opto per each rep
    flash_off_opto_long_vals=sum(spikes.*flopoff_long_array,2);
    flash_off_opto_long_vals=flash_off_opto_long_vals(flopoff_long_dur_reps>0);
    flash_off_opto_long_vals=flash_off_opto_long_vals./flopoff_long_dur_reps(flopoff_long_dur_reps>0); %firing rate per sec
        
    %shuffle: flash on vs. flash on+long opto
    [pval_s, mean_diff, effsize]  = permutationTest_with_subsampling(flash_off_vals, flash_off_opto_long_vals,n_tests);
    params_vals(ci).flopoff_long_meanfr=mean(flash_off_opto_long_vals,"omitnan");
    params_vals(ci).flopoff_long_s_pval=pval_s;    
    params_vals(ci).flopoff_long_s_effsize=effsize; 
    params_vals(ci).flopoff_long_s_obsdiff=mean_diff;
end

save(params_file,'params_vals');
end



function [rowi, sti]=get_bl_samples(b_arr,nsamples,win_r)
repinds = find(sum(b_arr,2)>win_r);
nreps=length(repinds);
rowiind=randi(nreps,nsamples,1);
rowi=repinds(rowiind);
tstart=arrayfun(@(i) find(b_arr(repinds(i),:),1,"first"),1:nreps);
tend=arrayfun(@(i) find(b_arr(repinds(i),:),1,"last"),1:nreps);
tend=tend-win_r;
bldur=tend-tstart;
allsample_dur=bldur(rowi);
sti = arrayfun(@(i) randi(allsample_dur(i)), 1:nsamples);
sti=sti+tstart(rowi)-1;
end


function fr_samples = get_fr_samples(spikes_arr, win_r, row_samples, row_t_index, sampling_rate)
nsamples=length(row_samples);
ent=row_t_index+win_r-1;
fr_samples = arrayfun(@(i) sum(spikes_arr(row_samples(i),row_t_index(i):ent(i))), 1:nsamples);
fr_samples = fr_samples*sampling_rate/win_r;
end
