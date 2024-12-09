function data  = set_csda_depth_per_cluster(data)
    %appr depth of the part of the Superficial Superior Culiculus above the
    %inflection point
    depth_SC_abov_inflection=300;

    %from cluster property file find the depth of the inflection channel
    if ~isfield(data,'csda')
        return;
    end
    if ~isfield(data,'spike_data')
        return;
    end

    %load channel map;
    load('chanMapCamNeu_Lin-A32-25_August2022_5.mat');
    cha_labels=flipud(chanMap);
    cha_labels=cha_labels-1;
    nchannels = numel(cha_labels);
    %create a map of channel labes and the SC depth determined thru CSDA
    cha_depth_SC=zeros(nchannels,2);
    cha_depth_SC(:,1)=cha_labels;
    cha_inflection = data.csda.inflection_cha_above;
    cha_inflection_id = find(cha_labels==cha_inflection);
    cha_spacing = data.csda.cha_spacing_um;
    for ch_i=cha_inflection_id:-1:1
        cha_depth_SC(ch_i,2)=depth_SC_abov_inflection-(cha_inflection_id-ch_i)*cha_spacing;
    end
    for ch_i=cha_inflection_id+1:nchannels
        cha_depth_SC(ch_i,2)=depth_SC_abov_inflection+(ch_i-cha_inflection_id)*cha_spacing;
    end
    %assign SC depth to the clusters
    nclusters = length(data.spike_data);
    for ci=1:nclusters
        cha_i=find(cha_depth_SC(:,1)==data.spike_data(ci).channel);
        data.spike_data(ci).SC_depth = cha_depth_SC(cha_i,2);
    end
end