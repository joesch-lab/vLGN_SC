function [forward_cm_s, sideways_cm_s, rotation_deg_s, ball_time] = ballcam2speed(ball, sample_rate)
% ballcam data is a 0/1 pulse signal
% the ball displacement is encoded by the length of the pulse

%% code below is an excefrom Toni's script 
%transform the length of 0-1 signal into the values [-127.127]

%from all ballcam get all onsets
onsets=[];
for b=1:4
%     onsets_temp=find([0 diff(ball{b})]==1);
%     offsets=find([0 diff(ball{b})]==-1);
    onsets_temp=find(diff([0 ball{b}])==1);
    offsets=find(diff([ball{b} 0])==-1);
%     if ball{b}(1)==1; offsets(1)=[];end
%     if ball{b}(end)==1; onsets_temp(end)=[];end
    onsets_temp=onsets_temp(1:numel(offsets));%/sample_rate;
    onsets=cat(2,onsets,onsets_temp);
end
onsets=unique(onsets); %get all onsets
onsets=onsets(diff(onsets)>2);

if isempty(onsets)
    forward_cm_s=[];
    sideways_cm_s=[];
    rotation_deg_s=[];
    ball_time = [];
    warning('No ballcam data found!');
    return;
end

onsets(1)=[];
%find periods where none of the channels have output 
% and add a missing onset
cycle_rate=median(diff(onsets));
d_onsets=diff(onsets);
%distance btw consequent onsets is larger than the cycle
onset_skips=find(d_onsets>cycle_rate*1.5);
onset_skip_num=round(d_onsets(onset_skips)/cycle_rate);
add_onsets=[];
% for all skips, insert missing trigger by interpolating btw neighbors
for n=1:numel(onset_skips)
    add_onsets=cat(2,add_onsets,round(linspace(onsets(onset_skips(n)),onsets(onset_skips(n)+1),onset_skip_num(n)+1)));
end
onsets=unique(cat(2,add_onsets,onsets));

%align all triggers, convert from pulse to pixels
ballC=cat(1,ball{:});
ball_data=nan(4,numel(onsets));
ball_time=onsets/sample_rate; %time in secs of all onsets
for n=1:numel(onsets)-1 %get value for ballcam per trigger
    %get the length of the ballcam signal per trigger
    temp=sum(ballC(:,onsets(n):onsets(n+1)-1),2);
    temp(temp>=127)=temp(temp>=127)-256;
    ball_data(:,n)=temp;
end
ball_data_raw=ball_data';
%remove nans
ball_time(any(isnan(ball_data_raw),2))=[];
ball_data_raw(any(isnan(ball_data_raw),2),:)=[];
ball_data=movmedian(ball_data_raw,5,1);% maybe only for calibration %smooth

%get weights to convert from pixels to cm
load('Calibration.mat');
Calibration=ModelParams.Calibration;

%transform from pixels to cm
forward_cm_s=sum(ball_data.*Calibration(1,:),2);
sideways_cm_s=sum(ball_data.*Calibration(2,:),2);
rotation_deg_s=sum(ball_data.*Calibration(3,:),2);

end