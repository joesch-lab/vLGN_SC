function transfer_CSDA_depth(datfile_nodepth, datfile_depth)
% transfers CSD information from the recording with the flash stimulus on which CSDA was performed
% (datfile_depth) to the recording of the same silicone probe position
% whithout CSD infomation (datfile_nodepth)
%
% O. Symonova 2024

dat1=load(datfile_nodepth);
if isfield(dat1.data.spike_data, 'SC_depth')
    disp([datfile_nodepth ' has already CSDA depth values']);
    return;
end
dat2=load(datfile_depth);
if ~isfield(dat2.data.spike_data, 'SC_depth')
    disp([datfile_depth ' has now CSDA depth values']);
    return;
end

%find the difference of depth between the cahnnel map and the CSDA
diff_depth=[];
for i=1:length(dat2.data.spike_data)
    diff_depth=[diff_depth, dat2.data.spike_data(i).depth-dat2.data.spike_data(i).SC_depth];
end
diff_depth=median(diff_depth); %it's always the same value

%apply the found difference to the recording without CSDA
for i=1:length(dat1.data.spike_data)
    dat1.data.spike_data(i).SC_depth=dat1.data.spike_data(i).depth-diff_depth;
end

data=dat1.data;
save(datfile_nodepth, 'data', '-v7.3');
end

%files to transfer depth information
% datfile_nodepth = 'D:\data\tnt_data\sin_checkers\Sham1-Gad2_recSC_M1F23-715_SinCheckers2_1768um_1\data.mat';
% datfile_depth = 'D:\data\tnt_data\Sham1-Gad2_recSC_M1F23-715_Olga_list_1768um_1_data.mat';

% datfile_nodepth = 'D:\data\tnt_data\sin_checkers\TNT-Gad2_MJ22-1842_Sin1_1\data.mat';
% datfile_depth = 'D:\data\tnt_data\TNT-Gad2_MJ22-1842_OlgaTestList_10_data.mat';

% datfile_nodepth = 'D:\data\tnt_data\sin_checkers\TNT-Gad2_MJ22-1842_Sin2_1\data.mat';
% datfile_depth = 'D:\data\tnt_data\TNT-Gad2_MJ22-1842_OlgaTestList_10_data.mat';