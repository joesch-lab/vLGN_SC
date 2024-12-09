function foldername_data = extract_data_from_foldername(fname)
    foldername_data={};
    foldername_data.folder = fname;
    [~,exp_name,~]=fileparts(fname);
    foldername_data.date=exp_name(1:6);

    return;
    
    %split the folder name with _ and ,
    strparts=strsplit(fname,{'_',','});

    %if there is a part 'withDi', Tomas used the dye    
    partid = find(contains(strparts,'withDi','IgnoreCase',true));
    if ~isempty(partid)
        foldername_data.withDye = 1;
    else
        foldername_data.withDye = 0;
    end
    
    %find the part that contains um, that the depth at the tip
    partid = find(contains(strparts,'um'));
    if isempty(partid)
        foldername_data.depth_at_tip = nan;
    else
        depthstr=strparts{partid};
        foldername_data.depth_at_tip = str2num(depthstr(1:end-2));
    end    

    %find the part that contains LPw, take the last recording
    partid = find(contains(strparts,'LPw','IgnoreCase',true));
     if isempty(partid)
        foldername_data.laser_power=nan;
    else
        foldername_data.laser_power=str2num(strparts{partid(end)+1});
    end   
    

    %find the part that contains Hz, take the last recording
    partid = find(contains(strparts,'Hz','IgnoreCase',true));
    if isempty(partid)
        foldername_data.opto_hz=nan;
    else
        hzstr=strparts{partid(end)};
        foldername_data.opto_hz=str2num(hzstr(1:end-2));
    end       
    

    %strarting from the end, collect all verbal parts - this is the note
    idnote=length(strparts);
    while ~isempty(regexpi(strparts{idnote},'[a-zA-Z]*'))    
        idnote=idnote-1;
    end
    if idnote<length(strparts)
        notst = strfind(fname,strparts(idnote+1));
        foldername_data.note=fname(notst:end);        
    else
        foldername_data.note='';
    end
end