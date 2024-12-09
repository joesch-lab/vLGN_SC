%% for all reps of sincheckers and all recordings, find and count saccades 
%motion of the pupil in azimuthal direction exceeding threshold. Compare control and TeLC groups
 
clear all; 
%Analysis OKR
load('okr_motion_data.mat');
still_tnt=[];
still_control=[];
move_tnt=[];
move_control=[];
rec_list={};
c=0;
c2=0;
saccade_threshold = 2.5;
for i=1:13        
    rec_list=cat(1,rec_list, okr_data(i).datfile);
    data_test = okr_data(i).paz;
    if okr_data(i).istnt ==1
        c=c+1;
        still = []; still2 =[]; mov =[]; mov2=[];
        still = data_test(mean(okr_data(i).movevec,2)==0,:);
        still2 = still-repmat(still(:,1), 1,41);
        mov = data_test(mean(okr_data(i).movevec,2)>0.5,:);
        mov2 = mov-repmat(mov(:,1), 1,41);
        still_tnt=[still_tnt; still2];
        move_tnt=[move_tnt; mov2];
        still_tnt_mean(c,:)=mean(still2);
        move_tnt_mean(c,:)=mean(mov2);

        %% calculate number of saccades in each session while still
        tnt_diff_temp =(diff(still2,2,2));
        tnt_diff_temp(tnt_diff_temp<saccade_threshold & tnt_diff_temp>-saccade_threshold )=0;
        tnt_diff_temp(tnt_diff_temp~=0)=1;
        temp=[];
        temp = reshape(tnt_diff_temp',1,[]);
        for i=1:length(temp)-1
            if temp(i) && temp(i+1)~=0
                temp(i+1) = 0;
            end
        end
        Saccade_count_tnt_still(c,:)=sum(temp)./size(still2,1);
        count(c,:)=size(still2,1);
        %% calculate number of saccades in each session shile moving
        tnt_diff_temp_mov =(diff(mov2,2,2));
        tnt_diff_temp_mov(tnt_diff_temp_mov<saccade_threshold & tnt_diff_temp_mov>-saccade_threshold)=0;
        tnt_diff_temp_mov(tnt_diff_temp_mov~=0)=1;
        temp=[];
        temp = reshape(tnt_diff_temp_mov',1,[]);
        for i=1:length(temp)-1
            if temp(i) && temp(i+1)~=0
                temp(i+1) = 0;
            end
        end
        Saccade_count_tnt_mov(c,:)=sum(temp)./size(mov2,1);

    else
        c2=c2+1;
        still = []; still2 =[]; mov =[]; mov2=[];
        still = data_test(mean(okr_data(1).movevec,2)==0,:);
        still2 = still-repmat(still(:,1), 1,41);
        mov = data_test(mean(okr_data(1).movevec,2)>0.1,:);
        mov2 = mov-repmat(mov(:,1), 1,41);
        still_control=[still_control; still2];
        move_control=[move_control; mov2];
        still_control_mean(c2,:)=mean(still2);
        move_control_mean(c2,:)=mean(mov2);

        %% calculate number of saccades in each session while still
        control_diff_temp =(diff(mov2,2,2));
        control_diff_temp(control_diff_temp<2 & control_diff_temp>-2 )=0;
        control_diff_temp(control_diff_temp~=0)=1;
        temp=[];
        temp = reshape(control_diff_temp',1,[]);
        for i=1:length(temp)-1
            if temp(i) && temp(i+1)~=0
                temp(i+1) = 0;
            end
        end
        Saccade_count_control_still(c2,:)=sum(temp)./size(still2,1);

        %% calculate number of saccades in each session shile moving
        control_diff_temp_mov =(diff(mov2,2,2));
        control_diff_temp_mov(control_diff_temp_mov<2 & control_diff_temp_mov>-2 )=0;
        control_diff_temp_mov(control_diff_temp_mov~=0)=1;
        temp=[];
        temp = reshape(control_diff_temp_mov',1,[]);
        for i=1:length(temp)-1
            if temp(i) && temp(i+1)~=0
                temp(i+1) = 0;
            end
        end
        Saccade_count_control_mov(c2,:)=sum(temp)./size(mov2,1);

    end
end


still_tnt_mean = mean(still_tnt(:,10:31),2);
[B, I]=sort(still_tnt_mean,'descend');
still_tnt = still_tnt(I,:);

still_control_mean = mean(still_control(:,10:31),2);
[B, I]=sort(still_control_mean,'descend');
still_control = still_control(I,:);


move_tnt_mean = mean(move_tnt(:,10:31),2);
[B, I]=sort(move_tnt_mean,'descend');
move_tnt = move_tnt(I,:);

move_control_mean = mean(move_control(:,10:31),2);
[B, I]=sort(move_control_mean,'descend');
move_control = move_control(I,:);


hist_still_tnt= hist(abs(diff(still_tnt(:))),0:0.1:10)./sum(hist(abs(diff(still_tnt(:))),0:0.1:10));
hist_still_control= hist(abs(diff(still_control(:))),0:0.1:10)./sum(hist(abs(diff(still_control(:))),0:0.1:10));

hist_move_tnt= hist(abs(diff(move_tnt(:))),0:0.1:10)./sum(hist(abs(diff(move_tnt(:))),0:0.1:10));
hist_move_control= hist(abs(diff(move_control(:))),0:0.1:10)./sum(hist(abs(diff(move_control(:))),0:0.1:10));

%Direction of saccade size
tnt_diff =(diff(still_tnt,2,2));
control_diff=(diff(still_control,2,2));
tnt_diff(tnt_diff<saccade_threshold & tnt_diff>-saccade_threshold )=0;
tnt_diff(tnt_diff>saccade_threshold)=1;
tnt_diff(tnt_diff<-saccade_threshold)=-1;
control_diff(control_diff<saccade_threshold & control_diff>-saccade_threshold )=0;
control_diff(control_diff>saccade_threshold)=1;
control_diff(control_diff<-saccade_threshold)=-1;

temp = reshape(control_diff',1,[]);
for i=1:length(temp)-1
    if temp(i) & temp(i+1)~=0
        temp(i+1) = 0;
    end
end
control_diff2 = reshape(temp',size(control_diff'));
control_diff =control_diff2';

temp = reshape(tnt_diff',1,[]);
for i=1:length(temp)-1
    if temp(i) & temp(i+1)~=0
        temp(i+1) = 0;
    end
end
tnt_diff2 = reshape(temp',size(tnt_diff'));
tnt_diff =tnt_diff2';

tnt_diff_p = tnt_diff;
tnt_diff_n = tnt_diff;
tnt_diff_p(tnt_diff_p<1)=0;
tnt_diff_n(tnt_diff_n>-1)=0;
tnt_diff_n=-tnt_diff_n;

mean_tnt_diff_p=mean(tnt_diff_p);
mean_tnt_diff_n=mean(tnt_diff_n);
mean_tnt_diff_pr = imresize(mean_tnt_diff_p,[1 8])./size(tnt_diff,1);
mean_tnt_diff_nr = imresize(mean_tnt_diff_n,[1 8])./size(tnt_diff,1);


control_diff_p = control_diff;
control_diff_n = control_diff;
control_diff_p(control_diff_p<1)=0;
control_diff_n(control_diff_n>-1)=0;
control_diff_n=-control_diff_n;
mean_control_diff_p=mean(control_diff_p);
mean_control_diff_n=mean(control_diff_n);
mean_control_diff_pr = imresize(mean_control_diff_p,[1 8])./size(control_diff,1);
mean_control_diff_nr = imresize(mean_control_diff_n,[1 8])./size(control_diff,1);

norm_value=max([mean_tnt_diff_nr mean_tnt_diff_nr mean_control_diff_pr mean_control_diff_nr]);
mean_control_diff_nr = mean_control_diff_nr./norm_value;
mean_control_diff_pr = mean_control_diff_pr./norm_value;
mean_tnt_diff_nr = mean_tnt_diff_nr./norm_value;
mean_tnt_diff_pr = mean_tnt_diff_pr./norm_value;

%% Direction of Saccade While Running
%Direction of saccade size
tnt_diffm =(diff(move_tnt,2,2));
control_diffm=(diff(move_control,2,2));
tnt_diffm(tnt_diffm<saccade_threshold & tnt_diffm>-saccade_threshold )=0;
tnt_diffm(tnt_diffm>saccade_threshold)=1;
tnt_diffm(tnt_diffm<-saccade_threshold)=-1;
control_diffm(control_diffm<saccade_threshold & control_diffm>-saccade_threshold )=0;
control_diffm(control_diffm>saccade_threshold)=1;
control_diffm(control_diffm<-saccade_threshold)=-1;

temp = reshape(control_diffm',1,[]);
for i=1:length(temp)-1
    if temp(i) & temp(i+1)~=0
        temp(i+1) = 0;
    end
end
control_diffm2 = reshape(temp',size(control_diffm'));
control_diffm =control_diffm2';

temp = reshape(tnt_diffm',1,[]);
for i=1:length(temp)-1
    if temp(i) & temp(i+1)~=0
        temp(i+1) = 0;
    end
end
tnt_diffm2 = reshape(temp',size(tnt_diffm'));
tnt_diffm =tnt_diffm2';

tnt_diff_pm = tnt_diffm;
tnt_diff_nm = tnt_diffm;
tnt_diff_pm(tnt_diff_pm<1)=0;
tnt_diff_nm(tnt_diff_nm>-1)=0;
tnt_diff_nm=-tnt_diff_nm;

mean_tnt_diff_pm=mean(tnt_diff_pm);
mean_tnt_diff_nm=mean(tnt_diff_nm);
mean_tnt_diff_prm = imresize(mean_tnt_diff_pm,[1 8])./size(tnt_diffm,1);
mean_tnt_diff_nrm = imresize(mean_tnt_diff_nm,[1 8])./size(tnt_diffm,1);


control_diff_pm = control_diffm;
control_diff_nm = control_diffm;
control_diff_pm(control_diff_pm<1)=0;
control_diff_nm(control_diff_nm>-1)=0;
control_diff_nm=-control_diff_nm;
mean_control_diff_pm=mean(control_diff_pm);
mean_control_diff_nm=mean(control_diff_nm);
mean_control_diff_prm = imresize(mean_control_diff_pm,[1 8])./size(control_diffm,1);
mean_control_diff_nrm = imresize(mean_control_diff_nm,[1 8])./size(control_diffm,1);

norm_valuem=max([mean_tnt_diff_nrm mean_tnt_diff_nrm mean_control_diff_prm mean_control_diff_nrm]);
mean_control_diff_nrm = mean_control_diff_nrm./norm_valuem;
mean_control_diff_prm = mean_control_diff_prm./norm_valuem;
mean_tnt_diff_nrm = mean_tnt_diff_nrm./norm_valuem;
mean_tnt_diff_prm = mean_tnt_diff_prm./norm_valuem;

%% Plot Heatmap of Saccades - Not running

figure
subplot(1,2,1)
imagesc(-control_diff)
title('Control OKR Saccades Still');
subplot(1,2,2)
imagesc(-tnt_diff)
title('TeLC OKR Saccades Still');
clim([-1 1]);
colormap(flipud(othercolor('RdBu7')));



%% plot heat map saccades still
figure
control_diffm(control_diffm~=0)=1;
tnt_diffm(tnt_diffm~=0)=1;
subplot(1,2,1)
imagesc(-control_diffm)
title('Control OKR Saccades move');
subplot(1,2,2)
imagesc(-tnt_diffm)
title('TeLC OKR Saccades move');
colormap('gray');

%% plot heat map saccades still
figure
control_diff(control_diff~=0)=1;
tnt_diff(tnt_diff~=0)=1;
subplot(1,2,1)
imagesc(-control_diff)
title('Control OKR Saccades still');
subplot(1,2,2)
imagesc(-tnt_diff)
title('TeLC OKR Saccades still');
colormap('gray');


%% plot heatmap of Eye motions - not running
figure
subplot(3,2,1)
imagesc(still_control)
title('control OKR Still')
clim([-5 5]);
colormap(flipud(othercolor('RdBu7')));
subplot(3,2,2)
imagesc(still_tnt)
title('TeLC OKR Still');
clim([-5 5]);
subplot(3,2,3)
plot((cos(pi:2*pi/41:3*pi)+1)./2,'b')
hold on
plot(mean(still_control),'k')
axis([1 41 0 1.5]);
subplot(3,2,4)
plot((cos(pi:2*pi/41:3*pi)+1)./2,'b')
hold  on
plot(mean(still_tnt),'r')
axis([1 41 0 1.5]);
subplot(3,2,5)
plot(mean_control_diff_pr,'k')
hold on
plot(mean_control_diff_nr,'k:')
axis([1 8 0 1.2]);
subplot(3,2,6)
plot(mean_tnt_diff_pr,'k')
hold on
plot(mean_tnt_diff_nr,'k:')
axis([1 8 0 1.2]);

mean_still_control=mean(still_control);
mean_still_tnt=mean(still_tnt);
stim_curve=(cos(pi:2*pi/40:3*pi)+1)./2;
timevec=linspace(0,4,41);
control_saccades_still=control_diff~=0;
tnt_saccades_still=tnt_diff~=0;
mean_control_sacc_nasal=mean_control_diff_pr;
mean_control_sacc_temporal=mean_control_diff_nr;
mean_tnt_sacc_nasal=mean_tnt_diff_pr;
mean_tnt_sacc_temporal=mean_tnt_diff_nr;
Saccade_count_tnt_still_mean=mean(Saccade_count_tnt_still);
Saccade_count_control_still_mean=mean(Saccade_count_control_still);
ranksumtest={};
[ranksumtest.p,ranksumtest.h,ranksumtest.stats]=ranksum(Saccade_count_tnt_still, Saccade_count_control_still);

 save('saccades_OKR.mat','still_control','still_tnt', 'mean_still_control','mean_still_tnt','stim_curve','timevec',...
'control_saccades_still','tnt_saccades_still','mean_control_sacc_nasal','mean_control_sacc_temporal','mean_tnt_sacc_nasal','mean_tnt_sacc_temporal',...
'Saccade_count_tnt_still','Saccade_count_control_still','Saccade_count_tnt_still_mean', 'Saccade_count_control_still_mean','ranksumtest','rec_list');

%% plot average still
figure
h=plot(ones(length(Saccade_count_tnt_still)), Saccade_count_tnt_still,'ro')
set(h, {'MarkerFaceColor'}, get(h,'Color')); 
hold on 
bar(1, mean(Saccade_count_tnt_still),'FaceColor','none', 'EdgeColor','r');
hold on 
h2=plot(2.*ones(length(Saccade_count_control_still)), Saccade_count_control_still,'ko')
set(h2, {'MarkerFaceColor'}, get(h2,'Color')); 
bar(2, mean(Saccade_count_control_still),'FaceColor','none', 'EdgeColor','k');
title(['ranksumtest pvalue=',num2str(ranksumtest.p)])
