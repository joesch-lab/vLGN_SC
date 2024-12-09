%% run this script to collect preprocessed RFs, 
%% parametrize them and compare RF with/without opto
% only stimuli with opto pulse <=100ms are selected (avoid pupil dilation)
% first collect signal-to-noise (s2n) variance from all units
% set s2n threshold at 80% percentile to select high-quality RFs
% parametrize RFs with the s2n> threshold
% for the the selected RFs, collect width and amplitude for opto +/-
% conditions
% make summary plot



% for all RF files find an appropriate snr value
resfolder = 'D:\data\resub\data\vertical_bars_opto';
dirlist=dir(fullfile(resfolder,'*rfinfo.mat'));
frame_history=30;

% use dates to find animal id for opto bar recordings
rf_animal_dates(1).id='MJ21_510';
rf_animal_dates(1).st_rec=datevec('210713','yymmdd');
rf_animal_dates(1).en_rec=datevec('210729','yymmdd');
rf_animal_dates(2).id='MJ21_508';
rf_animal_dates(2).st_rec=datevec('210903','yymmdd');
rf_animal_dates(2).en_rec=datevec('210910','yymmdd');
rf_animal_dates(3).id='MJ21_762';
rf_animal_dates(3).st_rec=datevec('210916','yymmdd');
rf_animal_dates(3).en_rec=datevec('210917','yymmdd');
rf_animal_dates(4).id='MJ21_763';
rf_animal_dates(4).st_rec=datevec('210923','yymmdd');
rf_animal_dates(4).en_rec=datevec('210924','yymmdd');

%% run the loop to collect all snr to build a distribution and select a threshold
s2n_all=[];
% check the distribution of snr and values for the recordings with t_opto <=0.1
for di=1:length(dirlist)
    if contains(dirlist(di).name,'no_laser', 'IgnoreCase',true), continue, end;

    rfresfile=fullfile(resfolder, dirlist(di).name);
    load(rfresfile);
    if stimparam.period_on>0.1
        continue;
    end

    for ci=1:length(RF)
        if isnan(RF(ci).snr_var_opto), continue; end       
        s2n_all=[s2n_all,RF(ci).snr_var_opto, RF(ci).snr_var_noopto];       
    end
end
snr_threhold  = quantile(s2n_all,0.8);

% [hc, be]=histcounts(s2n_all,100);
% bc=be(1:end-1)+diff(be(1:2))/2;
% fig=figure; bar(bc, hc);
% title('S2N for all RF');


%% for clusters with at least one rf of snr> threshold: collect snr +/-opto
%% for clusters with both rf of snr> threshold: fit gaussian
rffigname=fullfile(resfolder, 'randombar_rf_sel.pdf');
rf_bars={};
k=1;
gaussian_1d = @(x,data)x(1)*exp(-(((data-x(2)).^2)/(2*x(3).^2)))+x(4);

for di=1:length(dirlist)
    if contains(dirlist(di).name,'no_laser', 'IgnoreCase',true), continue, end;

    rfresfile=fullfile(resfolder, dirlist(di).name);
    load(rfresfile);
    if stimparam.period_on>0.1
        continue;
    end

    ncl=length(RF);

    for ci=1:ncl
        if isnan(RF(ci).logp2n_opto), continue; end

        if RF(ci).snr_var_noopto<snr_threhold && RF(ci).snr_var_opto<snr_threhold
            continue;
        end

        maxval=max(max(RF(ci).RFopto(:)),max(RF(ci).RFnoopto(:)));
        %time of max variance
        var_t=var(RF(ci).RFopto(1:frame_history,x1:x2),[],2);
        [~,t_maxvar_opto]=max(var_t);
        var_t=var(RF(ci).RFnoopto(1:frame_history,x1:x2),[],2);
        [~,t_maxvar_noopto]=max(var_t);

        rf_bars(k).opto.snr=RF(ci).snr_var_opto;
        rf_bars(k).noopto.snr=RF(ci).snr_var_noopto;
        rf_bars(k).opto.peak2noiselog=RF(ci).logp2n_opto;
        rf_bars(k).noopto.peak2noiselog=RF(ci).logp2n_noopto;
        rf_bars(k).spikes_used=RF(ci).numspikes_noopto_used;
        rf_bars(k).recname=dirlist(di).name;
        rf_bars(k).cl_id=ci;

        if RF(ci).snr_var_opto>=snr_threhold
            [g0, gfit, resnorm, rf_profile, rf_dynamics] = parametrize_RF(RF(ci).RFopto(:,x1:x2), frame_history, t_maxvar_opto);
            rf_bars(k).opto.t_maxvar=t_maxvar_opto;
            rf_bars(k).opto.data_estimate=g0;
            rf_bars(k).opto.fit=gfit;
            rf_bars(k).opto.fitresnorm=resnorm;
            rf_bars(k).opto.rf_profile=rf_profile;
            rf_bars(k).opto.rf_dynamics=rf_dynamics;
            rf_bars(k).opto.t_maxvar=t_maxvar_opto;
        end
        if RF(ci).snr_var_noopto>=snr_threhold
            [g0, gfit, resnorm, rf_profile, rf_dynamics] = parametrize_RF(RF(ci).RFnoopto(:,x1:x2), frame_history, t_maxvar_noopto);
            rf_bars(k).noopto.t_maxvar=t_maxvar_opto;
            rf_bars(k).noopto.data_estimate=g0;
            rf_bars(k).noopto.fit=gfit;
            rf_bars(k).noopto.fitresnorm=resnorm;
            rf_bars(k).noopto.rf_profile=rf_profile;
            rf_bars(k).noopto.rf_dynamics=rf_dynamics;
            rf_bars(k).noopto.t_maxvar=t_maxvar_noopto;
        end

        %make figure for fits of rf where comparison makes sense
        if RF(ci).snr_var_noopto>=snr_threhold && RF(ci).snr_var_opto>=snr_threhold

            fig=figure('Position',[100,100,800,600]);
            xspace=1:numel(rf_bars(k).opto.rf_profile);
            t=tiledlayout(5,2,"TileSpacing","compact");

            nexttile([3,1]);
            imagesc(RF(ci).RFopto(:,x1:x2), [0,maxval]);
            colormap(parula);
            hold on;
            plot([0.5,size(RF(ci).RFopto,2)+0.5],[frame_history+0.5 frame_history+0.5],'--w');
            title({'+ opto'; ['p2n=',num2str(RF(ci).logp2n_opto)];...
                ['snrvar=',num2str(RF(ci).snr_var_opto)]});

            nexttile([3,1]);
            imagesc(RF(ci).RFnoopto(:,x1:x2), [0,maxval]);
            colormap(parula);
            hold on;
            plot([0.5,size(RF(ci).RFnoopto,2)+0.5],[frame_history+0.5 frame_history+0.5],'--w');
            title({'- opto'; ['p2n=',num2str(RF(ci).logp2n_noopto)];...
                ['snrvar=',num2str(RF(ci).snr_var_noopto)]});

            %opto RF
            fitgaus=gaussian_1d(rf_bars(k).opto.fit,xspace);
            nexttile; plot(rf_bars(k).opto.rf_profile, 'k'); hold on;
            plot(fitgaus,'b');
            ylabel('spatial RFs');
            xticks([]);

            %no opto RF
            fitgaus=gaussian_1d(rf_bars(k).noopto.fit,xspace);
            nexttile; plot(rf_bars(k).noopto.rf_profile, 'k'); hold on;
            plot(fitgaus,'b');
            xticks([]);

            %temporal RF
            nexttile([1,2]);
            h1=plot(rf_bars(k).noopto.rf_dynamics, 'k'); hold on;
            h2=plot(rf_bars(k).opto.rf_dynamics, 'm');
            minval=min(min(rf_bars(k).noopto.rf_dynamics),min(rf_bars(k).opto.rf_dynamics));
            maxval=max(max(rf_bars(k).noopto.rf_dynamics),max(rf_bars(k).opto.rf_dynamics));
            plot([frame_history, frame_history],[minval,maxval],'--k');
            tmax=numel(rf_bars(k).noopto.rf_dynamics);
            xticks([1,frame_history-5, frame_history, tmax]);
            xticklabels(round([-frame_history/60, -5/60, 0, 5/60],2));
            xlabel('secs');
            ylabel('temporal RFs');
            legend([h1,h2],{'no opto', 'opto'},"Location","northwest");
            title(t,{rf_bars(k).recname; ['cluster ', num2str(rf_bars(k).cl_id)]; ['numspikes = ',num2str(RF(ci).numspikes_opto_used)]}, 'Interpreter', 'none');
            exportgraphics(fig, rffigname, 'BackgroundColor','none','ContentType','vector','Append', true);
            close(fig);
        end
        k=k+1;
    end
end
save(fullfile(resfolder, 'rf_bars_seledcted.mat'), 'rf_bars', 'snr_threhold',...
    's2n_opto', 's2n_noopto','-v7.3');


%% collect summary statistics: width and amplitude of RFs
figsummary=fullfile(resfolder, 'randombar_rf_summary.pdf');
w_opto=[];
w_noopto=[];
a_opto=[];
a_noopto=[];
reclist={};
j=1;
for ii=1:length(rf_bars)
    if rf_bars(ii).opto.snr>=snr_threhold && rf_bars(ii).noopto.snr>=snr_threhold
        reclist{j}=rf_bars(ii).recname;
        j=j+1;
        w_opto=[w_opto, rf_bars(ii).opto.fit(3)];
        w_noopto=[w_noopto,rf_bars(ii).noopto.fit(3)];
        a_opto=[a_opto, rf_bars(ii).opto.fit(1)];
        a_noopto=[a_noopto,rf_bars(ii).noopto.fit(1)];
    end
end
reclist=unique(reclist)';
anilist={};
j=1;
%get number of animals from the list of recordings
for ii=1:length(reclist)
    daterec=reclist{ii};
    daterec=datevec(daterec(1:6),'yymmdd');
    ki =  find(datetime(daterec)-datetime(cell2mat({rf_animal_dates(:).st_rec}'))>=0 &...
    datetime(daterec)-datetime(cell2mat({rf_animal_dates(:).en_rec}'))<=0);
    anilist{j}=rf_animal_dates(ki).id;
    j=j+1;
end
anilist=unique(anilist)';
% [amp_test.p, ~, amp_test.statistics] = signrank(a_noopto, a_opto);
% amp_test.testname='signrank';
%
% % maxa=max(max(a_opto),max(a_noopto));
% % fig= figure;
% % scatter(a_noopto, a_opto, 50, 'filled','MarkerFaceAlpha',0.3);
% % xlim([0,maxa]);
% % ylim([0,maxa]);
% % xlabel('no opto');
% % ylabel('opto');
% % hold on;
% % plot([0,maxa],[0,maxa],'--k');
% fig = figure; [sha,bha,a1] = scatterDiagHist(a_noopto, a_opto);
% sha.Marker='o';
% sha.MarkerFaceColor='flat';
% sha.MarkerEdgeAlpha=0.3;
% sha.MarkerFaceAlpha=0.3;
% bha.FaceAlpha=0.3;
% xlabel('amplitude - no opto');
% ylabel('amplitude  - opto');
% title({'strength of RF'; ['singrank test pvalue=',num2str(amp_test.p)];...
%     ['nn=', num2str(numel(a_opto)), '; nrec=', num2str(length(reclist)), '; nani=', num2str(length(anilist))]}, 'FontSize', 8);
% exportgraphics(fig,figsummary, 'BackgroundColor','none','ContentType','vector','Append', true);
%
pix2deg_scale=220/608; %horizontal span of the screen in deg and in pixels
[width_test.p, ~, width_test.statistics] = signrank(w_noopto*pix2deg_scale, w_opto*pix2deg_scale);
width_test.testname='signrank';
% maxw=max(max(w_opto),max(w_noopto));
% fig= figure;
% scatter(w_noopto, w_opto, 50, 'filled','MarkerFaceAlpha',0.3);
% xlim([0,maxw]);
% ylim([0,maxw]);
% xlabel('no opto');
% ylabel('opto');
% hold on;
% plot([0,maxw],[0,maxw],'--k');

fig = figure; [sha,bha,a1] = scatterDiagHist(w_noopto*pix2deg_scale, w_opto*pix2deg_scale);
sha.Marker='o';
sha.MarkerFaceColor='flat';
sha.MarkerEdgeAlpha=0.3;
sha.MarkerFaceAlpha=0.3;
bha.FaceAlpha=0.3;
xlabel('width - no opto (deg)');
ylabel('width  - opto (deg)');
title({'width of RF'; ['singrank test pvalue=',num2str(width_test.p)];...
    ['nn=', num2str(numel(w_opto)), '; nrec=', num2str(length(reclist)), '; nani=', num2str(length(anilist))]}, 'FontSize', 8);
exportgraphics(fig,figsummary, 'BackgroundColor','none','ContentType','vector','Append', true);


save(fullfile(resfolder, 'rf_bars_seledcted.mat'), 'rf_bars', 'snr_threhold',...
    's2n_opto', 's2n_noopto',...
    'reclist', 'anilist',...
    'a_noopto', 'a_opto', 'amp_test', ...
    'w_noopto', 'w_opto', 'width_test', 'pix2deg_scale', '-v7.3');



function [g0, gi, resnorm, rf_profile, rf_dynamics] = parametrize_RF(RFi, frame_history, maxvarT)
rf_win=1;
rfwidth=size(RFi,2);
RF_xspace=1:rfwidth;


%% parameters for the gaussian fit
myopts = optimoptions(@lsqcurvefit,'display','off','UseParallel',false);
gaussian_1d = @(x,data)x(1)*exp(-(((data-x(2)).^2)/(2*x(3).^2)))+x(4); %'a*exp(-((x-b)^2/(2c^2))+d';

tmin=max(1,maxvarT-rf_win);
tmax=min(frame_history,maxvarT+rf_win);
rf_profile=mean(RFi(tmin:tmax,:));

%fit Gaussian to the RF profile
[a0,c0]=max(abs(rf_profile));
sign_peak=sign(rf_profile(c0));
a0=a0*sign_peak;

w1=find((rf_profile(c0:-1:1)*sign_peak)<sign_peak*a0/2,1,'first');
if isempty(w1), w1=1; else, w1=c0-w1+1; end
w2=find((rf_profile(c0:end)*sign_peak)<sign_peak*a0/2,1,'first');
if isempty(w2), w2=rfwidth; else, w2=c0+w2-1; end
w_h2=w2-w1;
bl=median(rf_profile);

g0 =  [a0, c0, w_h2/2, bl];
[gi,resnorm,~,~] = lsqcurvefit(gaussian_1d ,g0, RF_xspace, rf_profile, [-Inf,1,1,-Inf],[Inf,rfwidth,rfwidth, Inf], myopts);
rf_dynamics=RFi(:,round(gi(2)));
end