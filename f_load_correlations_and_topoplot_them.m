clc
close all
warning('off','all')
addpath(genpath('/imaging/hp02/spm12b'));
addpath(genpath('/imaging/local/software/mne'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/new_functions');
%addpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/plotting');
addpath(genpath('/imaging/hp02/software_n_scripts/eeglab_current'));

%% load correlation data
%load('eeg/workspace_eeg_290616.mat');
%m = 3; % EEG
%load('megmag/workspace_megmag_290616.mat');
%m=1; %MEGMAG

%load('megcomb/workspace_161130_withlpf_with_non_para_stats.mat');
m=2; % MEGCOMB

% Individual participants:
%sname = [19 23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];

%% p threshold


thres_list = {'p05', 'p001', 'bonf'}; % Threshold options p05, p001, bonf
trialname = {'Action vs Scram', 'Spec vs Free'}; % 1 = spec, 2 = free
cor = {'all', 'acc'};
% Bonferroni correction:
%bonf = 0.05/(size(rvalues,1)*size(rvalues, 2)*2);
 

for thr = 1:3
    thres = thres_list{thr}
    for alloracc = 1:2
        allacc = cor{alloracc}
        clear thres_spec_free 
       
        %close all
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        
        
        base = '/imaging/hp02/finger_tapping08/analysis_spm/preprocessing/new_all_subjs/stimulus/ddm/lpf';
        
        prefix1 = 'Prm_acbefffMs'; % prefix of subject preprocessed files
     
        subj_file = sprintf('%s/s19/%s19_FT_blk1.mat', base, prefix1); % Subject file name

        cnt = 1;
        

        % load in SPM format to get sensor locations:
        D = spm_eeg_load(subj_file);
        modalities = {'MEGMAG', 'MEGMAG', 'EEG'}; % the second is really MEGCOMB, but the locations are not working for MEGCOMB for some reason
        labs = {'fTesla','fTesla/mm','\muV'};
        
        chans = find(strcmp(D.chantype,modalities{m}));
        
        
        cIndex = D.coor2D(D.indchantype(modalities{m}));
        Clabels = D.chanlabels(D.indchantype(modalities{m}));
        
        %remove bad channels
        if m == 3
            badchansEEGidx = find(ismember(chans,D.badchannels));
            goodchansEEG = setdiff(chans,D.badchannels);
            data = D(goodchansEEG,:,:);
            
            chLab = D.sensors(modalities{m}).label;
            chPos = D.sensors(modalities{m}).chanpos;
            
            chLab(badchansEEGidx)=[];
            chPos(badchansEEGidx,:)=[];
            axlim = [-8 8];
        else
            data = D(chans,:,:);
            sens = D.sensors('MEG');
            chLab = Clabels';
            [a,b,ind] = intersect(D.chanlabels(chans),sens.label);
            
            if m ==1 || m ==2
                chPos=sens.chanpos(ind,:);
                axlim = [-100 100];
            else
                %combine gradiometers
                rmsdata = reshape(data,2,102,size(data,2),size(data,3));
                rmsdata = squeeze(mean(rmsdata,1));
                %rmsdata(:,:,mu.badtrials)=[];
                data = rmsdata;
                chPos = sens.chanpos(ind(1:2:204),:);
                chLab = chLab(1:2:204);
                cIndex = cIndex(:,1:2:204);
                Clabels = Clabels(1:2:204);
                axlim = [-20 20];
                
            end
        end
        
        
        %% prepare for topo
        
       cd ./60-90Hz 
        
        meeglocs = [num2cell(1:length(chLab))',num2cell(chPos),num2cell(chLab)];
        
        dlmcell('plotlocations.xyz',meeglocs)
        clear EEG; EEG.nbchan = size(data,1);
        EEG.chanlocs = pop_chanedit(EEG, 'load',{'plotlocations.xyz' 'filetype' 'xyz'});
        
        % define data_vector
        data_vector = zeros(102,2);
        
        for ch  = 1:length(EEG.chanlocs)
            if alloracc == 1
                switch thres
                    case 'p05'
                        load('pvalues_p05_amp.mat')
                        
                        if pttest.all_action_scram(1,ch) ==1
                            data_vector(ch,1) = pttest.all_action_scram(3,ch);
                        end
                        if pttest.all_spec_free(1,ch) ==1
                            data_vector(ch,2) = pttest.all_spec_free(3,ch);
                        end
                        pthres(ch,1) = pttest.all_action_scram(1,ch);
                        pthres(ch,2) = pttest.all_spec_free(1,ch);
                        
                    case 'p001'
                        load('pvalues_p001_amp.mat')
                        
                        if pttest.all_action_scram(1,ch) ==1
                            data_vector(ch,1) = pttest.all_action_scram(3,ch);
                        end
                        if pttest.all_spec_free(1,ch) ==1
                            data_vector(ch,2) = pttest.all_spec_free(3,ch);
                        end
                        pthres(ch,1) = pttest.all_action_scram(1,ch);
                        pthres(ch,2) = pttest.all_spec_free(1,ch);
                    case 'bonf'
                        load('pvalues_pbonf_amp.mat')
                        
                        if pttest.all_action_scram(1,ch) ==1
                        data_vector(ch,1) = pttest.all_action_scram(3,ch);
                        end
                        if pttest.all_spec_free(1,ch) == 1
                            data_vector(ch,2) = pttest.all_spec_free(3,ch);
                        end
                        pthres(ch,1) = pttest.all_action_scram(1,ch);
                        pthres(ch,2) = pttest.all_spec_free(1,ch);
                        
                end
            else
                switch thres
                    case 'p05'
                        load('pvalues_p05_amp.mat')
                        if pttest.acc_action_scram(1,ch) ==1
                        data_vector(ch,1) = pttest.acc_action_scram(3,ch);
                        end
                        if pttest.acc_spec_free(1,ch) == 1
                            data_vector(ch,2) = pttest.acc_spec_free(3,ch);
                        end
                        pthres(ch,1) = pttest.acc_action_scram(1,ch);
                        pthres(ch,2) = pttest.acc_spec_free(1,ch);
                        
                    case 'p001'
                        load('pvalues_p001_amp.mat')
                        if pttest.acc_action_scram(1,ch) ==1
                        data_vector(ch,1) = pttest.acc_action_scram(3,ch);
                        end
                        if pttest.acc_spec_free(1,ch) == 1
                            data_vector(ch,2) = pttest.acc_spec_free(3,ch);
                        end
                        pthres(ch,1) = pttest.acc_action_scram(1,ch);
                        pthres(ch,2) = pttest.acc_spec_free(1,ch);
                        
                    case 'bonf'
                        load('pvalues_pbonf_amp.mat')
                        if pttest.acc_action_scram(1,ch) ==1
                        data_vector(ch,1) = pttest.acc_action_scram(3,ch);
                        end
                        if pttest.acc_spec_free(1,ch) == 1
                            data_vector(ch,2) = pttest.acc_spec_free(3,ch);
                        end
                        pthres(ch,1) = pttest.acc_action_scram(1,ch);
                        pthres(ch,2) = pttest.acc_spec_free(1,ch);
                        
                end
            end
        end
        
        titles = { 'Action vs Scrambled'; 'Specified vs Choice'; };
        
        for sp = 1:2
            subplot(1,2,sp)
            topoplot_hp(data_vector(:,sp), EEG.chanlocs, 'electrodes', 'off',...
                'style', 'fill');%, 'pmask', pthres(:,sp));
            title(titles(sp));
            caxis([-5 5]);
        end
        
        %cbar('vert',0,[0 1]);
        
         % save the figure
        figname = sprintf('%s_group_%s_nonpara_ttest_amp.fig', cor{alloracc}, thres_list{thr});
        set(gcf,'color', 'white');
        H = gcf;
        saveas(H, figname);
        close(H)
        clear H
        
        cd ..
    end
end
