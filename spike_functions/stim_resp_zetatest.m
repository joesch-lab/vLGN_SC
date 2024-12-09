function  neuron_zetavals = stim_resp_zetatest(tpts,spike_data, sr, stim_resp_dur)
%this is a wrapper for zeta-test for neurons in structure spike_data and
%stimulus event onset in array tpts
%spike timing can be stored seconds or in samples (will use sampling rate sr to convert)
neuron_zetavals={};

for ci=1:length(spike_data)
    spikes_all=[];
    if abs(sr-20000)<1
        spikest=double(spike_data(ci).spike_timing)/sr;
    else
        spikest=double(spike_data(ci).spike_timing);
    end
    try
        %% run the ZETA-test with default parameters
        %if we simply want to know if the neuron responds, no hassle, we can
        %use this simple syntax with default parameters:
        % dblUseMaxDur = min(diff(vecStimulusStartTimes)); %minimum of trial-to-trial durations
        intResampNum = 1000; %~50 random resamplings should give us a good enough idea if this cell is responsive, but if it's close to 0.05, we should increase this #. Generally speaking, more is better, so let's put 100 here.
        intPlot = 0;%what do we want to plot?(0=nothing, 1=inst. rate only, 2=traces only, 3=raster plot as well, 4=adds latencies in raster plot)
        vecRestrictRange = [0 inf];%do we want to restrict the peak detection to for example the time during stimulus? Then put [0 1] here.
        boolDirectQuantile = false;%if true; uses the empirical null distribution rather than the Gumbel approximation. Note that in this case the accuracy of your p-value is limited by the # of resamplings
        dblJitterSize = 1; %scalar value, sets the temporal jitter window relative to dblUseMaxDur (default: 2). Note that a value of "2" therefore means that the last data point that can be used is at t=last_onset + dblUseMaxDur*3
        boolStitch = false; %boolean, stitches data in the case of heterogeneous inter-event durations (default: true)

        %then run ZETA with those parameters
        [dblZetaP,sZETA,sRate,sLatencies] = zetatest(spikest,tpts,stim_resp_dur,intResampNum,intPlot,vecRestrictRange,boolDirectQuantile,dblJitterSize,boolStitch);

        %
        %     [dblZetaP,vecLatencies,sZETA,sRate] = getZeta(spikest, tpts, stim_resp_dur);
        neuron_zetavals(ci).z_pval=dblZetaP;
        neuron_zetavals(ci).z_frpeak_t = sRate.dblPeakTime;
        neuron_zetavals(ci).z_frpeak_maxfr = sRate.dblPeakRate;
        neuron_zetavals(ci).z_frpeak_w = sRate.dblPeakWidth;
        neuron_zetavals(ci).z_frpeak_tpts = sRate.vecT;
        neuron_zetavals(ci).z_ifr = sRate.vecRate;
    catch ME
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            disp('error: add zetatest to the matlab path');
        end
        neuron_zetavals(ci).z_pval=1;
        neuron_zetavals(ci).z_frpeak_t = nan;
        neuron_zetavals(ci).z_frpeak_maxfr = nan;
        neuron_zetavals(ci).z_frpeak_w =nan;
        neuron_zetavals(ci).z_frpeak_tpts = nan;
        neuron_zetavals(ci).z_ifr = nan;
    end
end

end