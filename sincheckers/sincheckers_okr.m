%% collect az and locomotion values for a specific number of repetitions
% use only one recording per day per animal. If multiple recordings per day
% are available use only the first

resfolder='D:\data\resub\data\sin_checkers';

%run sincheckers_analysis to create sindata in the file eye_vals_per_rep_all.mat
load(fullfile(resfolder, 'eye_vals_per_rep_all.mat'));

%how many reps to keep
selrepnum=225;
k=1;
processed_rec={}; %keep dates and animal ids of processed files
okr_data={};
%% collect az, move_yn vector, organize by rep, save if_tnt, if npx
%recordings in sindata are ordered by time, hence keep only the first
%recording per animal
for fi=1:length(sindata)     
    %skip if this is the second recording in the same day
    currdate=sindata(fi).date(1:11);
    if ~isempty(processed_rec)
        ii=find(ismember([processed_rec(:).ani_id],sindata(fi).ani_id));
        if ~isempty(ii) %multiple recordings per animal
            datesall=[processed_rec(ii).date];
            if contains(datesall,currdate) % same day recording
                continue;
            end
        end
    end   
    processed_rec(k).ani_id=sindata(fi).ani_id; %keep info of the proceed recording
    processed_rec(k).date=currdate;
    
    %time stamps of behavior frames
    minfrnum = min(numel(sindata(fi).camtrig_s), numel(sindata(fi).moveyn));
    camtrigs=sindata(fi).camtrig_s;

    %time vector of each rep
    repvals=sindata(fi).start_reps'+sindata(fi).bincen_t;
    
    %fit a sine wave to OKR
    mdl = fittype('a*sin(b*x+c)+d','indep','x');
    repcycle_dur_s=4;%sec
    lw= [0, 0.8*(2*pi)/repcycle_dur_s,-1.2*pi/2,-20];
    up= [40, 1.2*(2*pi)/repcycle_dur_s,-0.8*pi/2,20];
    xvals=linspace(0,4,41)'; 
    fitqual=zeros(selrepnum,1);
    fit_params=zeros(selrepnum,4);   
    %for each rep save locomotion classification values
    movevec=zeros(selrepnum,numel(sindata(fi).bincen_t));

    %for each rep fit sine function, collect locomotion data
    for ri=1:selrepnum
        yvals=sindata(fi).eye_vals_az(ri,:)';
        movevec(ri,:)=interp1(camtrigs(1:minfrnum),sindata(fi).moveyn(1:minfrnum),repvals(ri,:),'nearest');

        %fit sine function
        fit_mdl = fit(xvals,yvals,mdl,'start',[range(yvals)/2,(2*pi)/repcycle_dur_s,pi/2,mean(yvals)],'lower', lw,'upper',up);
        yfit = fit_mdl(xvals);
        SStot = sum((yvals-mean(yvals)).^2); % Total Sum-Of-Squares
        SSres = sum((yvals-yfit).^2);  % Residual Sum-Of-Squares
        qi=1-SSres/SStot;
        fitqual(ri)=qi;
        fit_params(ri,:)=[fit_mdl.a, fit_mdl.b, fit_mdl.c, fit_mdl.d];
    end
    okr_data(k).ani_id=sindata(fi).ani_id;
    okr_data(k).date=currdate;
    okr_data(k).recid=fi;
    okr_data(k).istnt=sindata(fi).istnt;
    okr_data(k).datfile=sindata(fi).datafile;
    okr_data(k).paz=sindata(fi).eye_vals_az(1:selrepnum,:);
    okr_data(k).movevec=movevec;
    okr_data(k).fit_params=fit_params;
    okr_data(k).fit_quality=fitqual;
    okr_data(k).is_npx=is_npx;
    k=k+1;
end
save(fullfile(resfolder, 'okr_motion_data.mat'),'okr_data', '-v7.3');

