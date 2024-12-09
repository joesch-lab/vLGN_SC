function [tstart, tend] = start_end_stimulus(redsignal, stimid, logfiles, force_intractive, sr)
%shows the plot of the red syncing signal (usually shown during the stimulation)
%if >1 stimulus, need to click on the empty space (no syncing signal)
%before and after each stimulus
%manual procedure because some intra-stimulus pauses are larger than
%intervals between stimuli

if ~exist("force_intractive","var") || isempty(force_intractive)
    force_intractive=0;
end

if ~exist("sr","var") || isempty(sr)
    sr=20000;
end

rft=find(redsignal);

if force_intractive==0 && length(logfiles)==1    
    tstart = rft(1);
    tend = rft(end);
    return;
end

padaround=120;

stimlogfile=fullfile(logfiles(stimid).folder, logfiles(stimid).name);
stimparam = read_stimulus_file(stimlogfile);
s1 = stimparam.tstart;
s2 = stimparam.tend;
% if ischar(stimparam.tstart)
%     s1 = datevec([stimparam.tstart],'yyyy/mm/dd HH:MM:SS:FFF');
%     s2 = datevec(stimparam.tend,'yyyy/mm/dd HH:MM:SS:FFF');
% else
%     s1 = datevec([stimparam.tstart],'yyyy/mm/dd HH:MM:SS:FFF');
%     s2 = datevec(stimparam.tend,'yyyy/mm/dd HH:MM:SS:FFF');
% end
stimdur = seconds(datetime(s2)-datetime(s1));

if stimid==1
    t0e=1;
    tNe=min(t0e+stimdur*sr+padaround*sr,length(redsignal));
elseif stimid==length(logfiles)
    tNe = length(redsignal);
    t0e =tNe - stimdur*sr-padaround*sr;
else
    stimlogfile1=fullfile(logfiles(1).folder, logfiles(1).name);
    stimparam1 = read_stimulus_file(stimlogfile1);
    s11 = stimparam1.tstart;
    s12 = stimparam1.tend;
%     s11 = datevec(stimparam1.tstart,'yyyy/mm/dd HH:MM:SS:FFF');
%     s12 = datevec(stimparam1.tend,'yyyy/mm/dd HH:MM:SS:FFF');
    tbtwstim = seconds(datetime(s1)-datetime(s11));
    tbtwstim_end = seconds(datetime(s2)-datetime(s11));
    t0e=round(max(1,tbtwstim*sr-padaround*sr));
    tNe=round(min(length(redsignal), tbtwstim_end*sr+padaround*sr));
end


disp('displaying the red frames around the stimulus');
disp(['Locate the gaps around the stimulus', logfiles(stimid).name]);
disp('list of stimuli');
for li=1:length(logfiles)
    disp(logfiles(li).name);
end
disp('now click on the GAPS around the stimulus and then press return');
scr_s = get(0,'screensize');
figure('Position',[50,50, scr_s(3)-100, 400]); plot(redsignal(uint32(t0e):uint32(tNe)));
[x,y]=ginput;
t0 = max(1,uint32(t0e+x(1)-1));
tN = min(length(redsignal),uint32(t0e+x(2)-1));

rft1=find(redsignal(t0:tN));
tstart=t0 + rft1(1)-1;
tend=t0 + rft1(end)-1;

figure, plot(redsignal);
hold on; scatter(tstart,1,'*r');
hold on; scatter(tend,1,'*r');

end

