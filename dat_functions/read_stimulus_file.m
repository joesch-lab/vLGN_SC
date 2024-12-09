function stimparams = read_stimulus_file(stimlogfile)
if endsWith(stimlogfile, '.mat')
    stimparams = load(stimlogfile);
elseif contains(stimlogfile, 'checker_stimuli_Jan_lowres')
    stimparams = read_checker_lowres(stimlogfile);
elseif contains(stimlogfile, 'checker_stimuli_Jan')
    stimparams = read_checker_stimuli_general_log(stimlogfile);
elseif contains(stimlogfile, 'stimuli_optomotor_gratings')
    stimparams = read_optomotor_gratings(stimlogfile);
elseif contains(stimlogfile, 'stimuli_full_field_curve_' + digitsPattern(6))
    stimparams = read_full_field_curve(stimlogfile);
elseif contains(stimlogfile, 'stimuli_full_field_' + digitsPattern(6))
    stimparams = read_full_field(stimlogfile);
else 
    disp(['No parser for ' stimlogfile]);
    stimparams=[];
    return;
end

if ischar(stimparams.tstart)
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';    
    stimparams.tstart=datevec(stimparams.tstart,formatIn);
    stimparams.tend=datevec(stimparams.tend,formatIn);
end
end


function stimparams = read_checker_lowres(stimlogfile)
    stimparams=[];
    if ~contains(stimlogfile,'checker_stimuli_Jan_lowres')
        return;
    end
        
    %% read the file with stimuli parameters
    f=fopen(stimlogfile);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

    % stimuli name line
    
    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';    
    stimparams.tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tend=datevec(tstr,formatIn);

    tline=fgetl(f);
    ifi_str=tline(length('Freq: ')+1:end);
    stimparams.ifi=str2double(ifi_str); %flip rate : frames per second

    tline=fgetl(f);
    t_str=tline(length('t\t')+1:end);
    stimparams.t=str2double(t_str); %requested time

    tline=fgetl(f);
    nfr_str=tline(length('nframes\t')+1:end);
    stimparams.nframes=str2num(nfr_str); %number of frames each random pattern was presented

    tline=fgetl(f);
    ncol_str=tline(length('ncol\t')+1:end);
    stimparams.ncol=str2num(ncol_str); %number of columns

    tline=fgetl(f);
    stimparams.distr=tline(length('distr\t')+1:end); %name of the distribution

    tline=fgetl(f);
    fraction_str=tline(length('fraction\t')+1:end);
    stimparams.fraction=str2double(fraction_str); %used in uniform distribution
    
    tline=fgetl(f);
    shift_fraction_str=tline(length('shift fraction\t')+1:end);
    stimparams.shift_fraction=str2double(shift_fraction_str); %for moving checkers    

    tline=fgetl(f);
    shift_steps_str=tline(length('shift steps\t')+1:end);
    stimparams.shift_steps=str2double(shift_steps_str); %steps for moving checkers    

    tline=fgetl(f);
    ri_tr=tline(length('maxRedIntensity\t')+1:end);
    stimparams.maxRedIntensity=str2num(ri_tr); %maximum red intencity

    tline=fgetl(f);
    ch_str=tline(length('channels\t[')+1:end-1);
    ch1=str2num(ch_str(1)); %used in uniform distribution
    ch2=str2num(ch_str(3));
    ch3=str2num(ch_str(5));
    stimparams.channels=[ch1 ch2 ch3];

    tline=fgetl(f);
    str_str=tline(length('stimulus frames presented: ')+1:end);
    stimparams.frames_presented=str2double(str_str); 
        
    tline=fgetl(f); %seed structure    
    tline=fgetl(f);
    seed_type=tline(length('Type\t')+1:end);
    tline=fgetl(f);
    seed_id= str2num(tline(length('Seed\t')+1:end));
    tline=fgetl(f);
    seed_size= str2num(tline(length('State (size)\t')+1:end));
    tline=fgetl(f);
    nn = strsplit(tline, {' '});
    nn= nn(~cellfun('isempty',nn));  
    nnum=uint32(cellfun(@str2num, nn))';    

    %random seed
    s.Type=seed_type;
    s.Seed = seed_id;
    s.State = nnum;
    stimparams.seed =s;

    fclose(f);
end

function stimparams = read_optomotor_gratings(stimlogfile)
    stimparams=[];
    if ~contains(stimlogfile,'optomotor_gratings')
        return;
    end
        
    %% read the file with stimuli parameters
    f=fopen(stimlogfile);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

    % stimuli name line
    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tend=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('Time per grating: ')+1:end);
    stimparams.t_grating=str2double(tstr);    

    tline=fgetl(f);
    ifi_str=tline(length('Flip rate: ')+1:end);
    stimparams.ifi=str2double(ifi_str);

    tline=fgetl(f);
    t_str=tline(length('Repetitions: ')+1:end);
    stimparams.nrep=str2double(t_str); 

    tline=fgetl(f);
    t_str=tline(length('Speed: ')+1:end);
    stimparams.speed_grating=str2double(t_str); 

    tline=fgetl(f);
    t_str=tline(length('With gray: ')+1:end);
    stimparams.with_gray=str2double(t_str); 
    
    tline=fgetl(f);
    t_str=tline(length('Spatial frequency: ')+1:end);
    stimparams.spfreq_grating=str2double(t_str);

    fclose(f);
end

function stimparams = read_full_field_curve(stimlogfile)
    stimparams=[];
    if ~contains(stimlogfile,'full_field_curve')
        return;
    end
        
    %% read the file with stimuli parameters
    f=fopen(stimlogfile);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

    % stimuli name line
    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tend=datevec(tstr,formatIn);

    tline=fgetl(f);
    ifi_str=tline(length('Flip rate: ')+1:end);
    stimparams.ifi=str2double(ifi_str);

    tline=fgetl(f);
    t_str=tline(length('Number of reps: ')+1:end);
    stimparams.nrep=str2double(t_str); 

    tline=fgetl(f);
    tstr=tline(length('tbefore: ')+1:end);
    stimparams.tbefore=str2double(tstr);    

    tline=fgetl(f);
    tstr=tline(length('tafter: ')+1:end);
    stimparams.tafter=str2double(tstr);  

    tline=fgetl(f);
    t_str=tline(length('Background_intensity: ')+1:end);
    stimparams.bgI =str2double(t_str); 

    tline=fgetl(f);
    ch_str=tline(length('Channels: [')+1:end-1);
    ch1=str2num(ch_str(1)); %used in uniform distribution
    ch2=str2num(ch_str(3));
    ch3=str2num(ch_str(5));
    stimparams.channels=[ch1 ch2 ch3];

    tline=fgetl(f);
    t_str=tline(length('Stimulus frames: ')+1:end);
    stimparams.nframes =str2double(t_str); 
    fclose(f);
end

function stimparams = read_full_field(stimlogfile)
    stimparams=[];
    if ~contains(stimlogfile,'full_field')
        return;
    end
        
    %% read the file with stimuli parameters
    f=fopen(stimlogfile);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

    % stimuli name line
    tline=fgetl(f);
    tstr=tline(length('Sart time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tend=datevec(tstr,formatIn);

    tline=fgetl(f);
    t_str=tline(length('Duration(sec): ')+1:end);
    stimparams.tdur=str2double(t_str);

    tline=fgetl(f);
    t_str=tline(length('Frames presented: ')+1:end);
    stimparams.nfr=str2double(t_str); 

    tline=fgetl(f);
    ifi_str=tline(length('Flip rate: ')+1:end);
    stimparams.ifi=str2double(ifi_str);   

    tline=fgetl(f);
    tstr=tline(length('Color: ')+1:end);
    stimparams.intensity=str2double(tstr);    

    tline=fgetl(f);
    ch_str=tline(length('Channel: [')+1:end-1);
    ch1=str2num(ch_str(1)); %used in uniform distribution
    ch2=str2num(ch_str(3));
    ch3=str2num(ch_str(5));
    stimparams.channels=[ch1 ch2 ch3];
    
    fclose(f);
end

function stimparams = read_checker_stimuli_general_log(stimlogfile)
    %% read the file with stimuli parameters
    f=fopen(stimlogfile);
    tline=fgetl(f);
    while  ~feof(f) && tline==""
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'checker'))
        disp('The file does not contain information about scanning dot stimulus');
        fclose(f);
        return;
    end

    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparams.tend=datevec(tstr,formatIn);

    tline=fgetl(f);
    ifi_str=tline(length('Freq: ')+1:end);
    stimparams.ifi=str2double(ifi_str); %flip rate : frames per second

    tline=fgetl(f);
    t_str=tline(length('t\t')+1:end);
    stimparams.t=str2double(t_str); %requested time

    tline=fgetl(f);
    res=strfind(tline, 'nframes');
    if isempty(res)
        nfr_str=tline(length('Frames presented:\t')+1:end);
        stimparams.frames_presented=str2num(nfr_str); %number of frames each random pattern was presented
        tline=fgetl(f);
    end

    nfr_str=tline(length('nframes\t')+1:end);
    stimparams.nframes=str2num(nfr_str); %number of frames each random pattern was presented

    tline=fgetl(f);
    ncol_str=tline(length('ncol\t')+1:end);
    stimparams.ncol=str2num(ncol_str); %number of columns

    tline=fgetl(f);
    stimparams.distr=tline(length('distr\t')+1:end); %name of the distribution

    tline=fgetl(f);
    fraction_str=tline(length('fraction\t')+1:end);
    stimparams.fraction=str2double(fraction_str); %used in uniform distribution
    
    tline=fgetl(f);
    shift_fraction_str=tline(length('shift fraction\t')+1:end);
    stimparams.shift_fraction=str2double(shift_fraction_str); %for moving checkers    

    tline=fgetl(f);
    ri_tr=tline(length('maxRedIntensity\t')+1:end);
    maxRedIntensity=str2num(ri_tr); %maximum red intencity

    tline=fgetl(f);
    ch_str=tline(length('channels\t[')+1:end-1);
    ch1=str2num(ch_str(1)); %used in uniform distribution
    ch2=str2num(ch_str(3));
    ch3=str2num(ch_str(5));
    stimparams.channels=[ch1 ch2 ch3];

    tline=fgetl(f);
    if ~isempty(strfind(tline,'stimulus frames'))
        str_str=tline(length('stimulus frames presented: ')+1:end);
        stimparams.frames_presented=str2double(str_str); 
        tline=fgetl(f);
    else
        stimparams.frames_presented=-1;
    end
    
    %% for the newer version with fadein and fadeout options
%     tline=fgetl(f); %this is either seed structure or fadein time
    fres=strfind(tline, 'Seed');
    if isempty(fres) %this is not seed structure, hence fadein
        fit_str=tline(length('fadein time: ')+1:end);
        stimparams.fadeintime=str2double(fit_str); %fadein time
        tline=fgetl(f); %fadeout time
        fot_str=tline(length('fadeout time: ')+1:end);
        stimparams.fadeouttime=str2double(fot_str); 
        tline=fgetl(f); % fadinframes
        fif_str=tline(length('fadein frames presented: ')+1:end);
        stimparams.fadeinframes=str2num(fif_str); 
        tline=fgetl(f); % stimulus frames
        stf_str=tline(length('stimulus frames presented: ')+1:end);
        stimparams.stimulusframes=str2num(stf_str);
        tline=fgetl(f); % fadeout frames
        fof_str=tline(length('fadeout frames presented: ')+1:end);
        stimparams.fadeoutframes=str2num(fof_str);    
    else    
        stimparams.fadeintime=0;   
        stimparams.fadeouttime=0; 
        stimparams.fadeinframes=0;     
        stimparams.stimulusframes=0;     
        stimparams.fadeoutframes=0; 
    end 

    %now read info about the random generator seed:
    if isempty(fres)
        tline=fgetl(f);
    end
    tline=fgetl(f);
    seed_type=tline(length('Type\t')+1:end);
    tline=fgetl(f);
    seed_id= str2num(tline(length('Seed\t')+1:end));
    tline=fgetl(f);
    seed_size= str2num(tline(length('State (size)\t')+1:end));
    tline=fgetl(f);
    nn = strsplit(tline, {' '});
    nn= nn(~cellfun('isempty',nn));  
    nnum=uint32(cellfun(@str2num, nn))';

    fclose(f);

    %restore the random seed
    s.Type=seed_type;
    s.Seed = seed_id;
    s.State = nnum;
    stimparams.s = s;
end