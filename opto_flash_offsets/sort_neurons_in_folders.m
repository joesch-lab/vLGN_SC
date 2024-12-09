
% load processed info about flash opto responses and sort the pictures
% into the folders to test the results of the tests

folder_root='D:\data\resub\data';
list_process = 'D:\data\resub\data\with_laser_list.txt'; %laser_yn=1;
folder_list = get_list_from_file(list_process);
path_to_figures= 'D:\data\resub\figs';
folder_group = 'D:\data\resub\figs\sort_neurons'

sSC_dir='sSC';
dSC_dir='dSC';
neuron_types={'vis_opto','vis_no_opto','opto_no_vis','rbnd','other'};
for ni=1:length(neuron_types)
    mkdir(fullfile(fullfile(folder_group,sSC_dir),neuron_types{ni}));
    mkdir(fullfile(fullfile(folder_group,dSC_dir),neuron_types{ni}));
end


% params for sSC border
minD=0;
maxD=Inf;
sSC_border=500;
%

reclist = get_list_from_file(list_process);
nf=length(reclist);
n_skipped=0;
pthresh=0.01;
for di=1:nf
    pref=reclist{di};
    data_dir=fullfile(folder_root,pref);

    disp(data_dir);
    ofdata_file = fullfile(data_dir,'opto_flash_data.mat');
    load(ofdata_file);
    param_file = fullfile(data_dir,'params_flash_opto_responses_all.mat');
    load(param_file);

    ncl=length(params_vals);
    for ci=1:ncl
        cl_im=fullfile(fullfile(data_dir,'cluster_im'),['cluster_',num2str(ci),'.png']);
        clim_new=[num2str(ci),'__',pref,'.png'];

        if of_data.cluster_depth(ci)>=minD
            if of_data.cluster_depth(ci)<sSC_border
                neuron_dir=fullfile(folder_group,sSC_dir);
            else
                neuron_dir=fullfile(folder_group,dSC_dir);
            end
        else
            n_skipped=n_skipped+1;
            continue; %depth <0 above SC
        end

        %% visual
        if is_visual(params_vals(ci),pthresh) %visual
            if is_opto(params_vals(ci),pthresh) % opto
                copyfile(cl_im,fullfile(fullfile(neuron_dir,'vis_opto'),clim_new));
                %check a rebound property
                if params_vals(ci).opoff_z_pval<=pthresh % rebound: time-locked event after opto-release
                    copyfile(cl_im,fullfile(fullfile(neuron_dir,'rbnd'),clim_new));
                end
            else %visual, no opto
                copyfile(cl_im,fullfile(fullfile(neuron_dir,'vis_no_opto'),clim_new));
            end
        else % non-visual
            if is_opto(params_vals(ci), pthresh) %non visual but opto
                copyfile(cl_im,fullfile(fullfile(neuron_dir,'opto_no_vis'),clim_new));
                %check a rebound property
                if params_vals(ci).opoff_z_pval<=pthresh % rebound: time-locked event after opto-release
                    copyfile(cl_im,fullfile(fullfile(neuron_dir,'rbnd'),clim_new));
                end
            else %no vis, no opto
                %non visual and non opto
                copyfile(cl_im,fullfile(fullfile(neuron_dir,'other'),clim_new));
            end
        end
    end
end
%% count the files in the sorted folders to count neurons
% sSC
sSC_fulldir=fullfile(folder_group,sSC_dir);
neuronN.sSC.n_other=length(dir(fullfile(fullfile(sSC_fulldir,'other'),'*.png')));
neuronN.sSC.n_vis_opto=length(dir(fullfile(fullfile(sSC_fulldir,'vis_opto'),'*.png')));
neuronN.sSC.n_vis_no_opto=length(dir(fullfile(fullfile(sSC_fulldir,'vis_no_opto'),'*.png')));
neuronN.sSC.n_opto_no_vis=length(dir(fullfile(fullfile(sSC_fulldir,'opto_no_vis'),'*.png')));
neuronN.sSC.n_rebound=length(dir(fullfile(fullfile(sSC_fulldir,'rbnd'),'*.png')));

neuronN.sSC.n_all=neuronN.sSC.n_vis_opto+neuronN.sSC.n_vis_no_opto+neuronN.sSC.n_opto_no_vis+neuronN.sSC.n_other;
neuronN.sSC.n_visual=neuronN.sSC.n_vis_opto+neuronN.sSC.n_vis_no_opto;
neuronN.sSC.n_opto=neuronN.sSC.n_vis_opto+neuronN.sSC.n_opto_no_vis;

% dSC
dSC_fulldir=fullfile(folder_group,dSC_dir);
neuronN.dSC.n_other=length(dir(fullfile(fullfile(dSC_fulldir,'other'),'*.png')));
neuronN.dSC.n_vis_opto=length(dir(fullfile(fullfile(dSC_fulldir,'vis_opto'),'*.png')));
neuronN.dSC.n_vis_no_opto=length(dir(fullfile(fullfile(dSC_fulldir,'vis_no_opto'),'*.png')));
neuronN.dSC.n_opto_no_vis=length(dir(fullfile(fullfile(dSC_fulldir,'opto_no_vis'),'*.png')));
neuronN.dSC.n_rebound=length(dir(fullfile(fullfile(dSC_fulldir,'rbnd'),'*.png')));

neuronN.dSC.n_all=neuronN.dSC.n_vis_opto+neuronN.dSC.n_vis_no_opto+neuronN.dSC.n_opto_no_vis+neuronN.dSC.n_other;
neuronN.dSC.n_visual=neuronN.dSC.n_vis_opto+neuronN.dSC.n_vis_no_opto;
neuronN.dSC.n_opto=neuronN.dSC.n_vis_opto+neuronN.dSC.n_opto_no_vis;

%all
neuronN.all.n_all=neuronN.sSC.n_all+neuronN.dSC.n_all;
neuronN.all.n_visual=neuronN.sSC.n_visual+neuronN.dSC.n_visual;
neuronN.all.n_opto=neuronN.sSC.n_opto+neuronN.dSC.n_opto;


neuronN.all.n_vis_opto=neuronN.sSC.n_vis_opto+neuronN.dSC.n_vis_opto;
neuronN.all.n_vis_no_opto=neuronN.sSC.n_vis_no_opto+neuronN.dSC.n_vis_no_opto;
neuronN.all.n_opto_no_vis=neuronN.sSC.n_opto_no_vis+neuronN.dSC.n_opto_no_vis;
neuronN.all.n_other=neuronN.sSC.n_other+neuronN.dSC.n_other;
neuronN.all.n_rebound=neuronN.sSC.n_rebound+neuronN.dSC.n_rebound;

neuronN.skipped = n_skipped;

filename=fullfile(path_to_figures,'neuron_numbers.mat');
save(filename,'neuronN','pthresh','reclist');



