clear all
clc
%warning off
%dbstop if error
%% Root OUTPUT DIRECTORY (e.g. '/imaging/meg.ryan/data')
bwd = '';

% Uncomment for the first run in MATLAB:
% addpath(genpath('/imaging/hp02/spm12b'));
% addpath(genpath('/imaging/local/software/mne'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/new_functions');

% File Containing CHANNEL INFORMATION for MEG and EEG (e.g. '/imaging/meg.ryan/Batch/chan_select_MEG_EEG_STI101.mat')
chanfile = '/imaging/hp02/mmn_08/analysis_spm/channels.mat';


% Specify (multiple) SUBJECTs)
cnt = 0;

%% Enter Subjects - Change Accordingly:
sname = [18 19 23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];
blkout = {'FT_blk1', 'FT_blk2', 'FT_blk3', 'FT_blk4'};
mf_path = '/imaging/hp02/finger_tapping08/maxFilt/';

subjects{1} = [mf_path, 'Lauras_mf/meg08_0018/MF'];
blocksin{1} = {'FT_trans_meg08_0018_1', 'FT_trans_meg08_0018_2', 'FT_trans_meg08_0018_3', 'FT_trans_meg08_0018_4'};   % as named after Maxfilter

for i = 2:13
    subjects{i} = sprintf('%sMF-rik/meg08_00%d', mf_path, sname(i));
    blocksin{i} = {'block1_raw_trans1stdef', 'block2_raw_trans1stdef', 'block3_raw_trans1stdef', 'block4_raw_trans1stdef'};   % as named after Maxfilter
    
end
% some exceptions:
blocksin{5} = {'block1_raw_trans1stdef', 'block2_raw_trans1stdef', 'block3_raw_trans1stdef', 'block_4_raw_trans1stdef'};
blocksin{7} = {'BLOCK_1_RAW_trans1stdef', 'BLOCK_2_RAW_trans1stdef', 'BLOCK3_RAW_trans1stdef', 'BLOCK_4_RAW_trans1stdef'};   % as named after Maxfilter
blocksin{8} = {'block1_raw_trans1stdef', 'block_2_raw_trans1stdef', 'block3_raw_trans1stdef', 'block_4_raw_trans1stdef'};   % as named after Maxfilter


for i = 14:19
    subjects{i} = sprintf('%sMF-rik/meg15_0%d', mf_path, sname(i));
    blocksin{i} = {'FT_blk1_raw_trans1stdef', 'FT_blk2_raw_trans1stdef', 'FT_blk3_raw_trans1stdef', 'FT_blk4_raw_trans1stdef'};   % as named after Maxfilter
end

% some exceptions:
blocksin{14} = {'blk1_raw_trans1stdef', 'blk2_raw_trans1stdef', 'blk3_raw_trans1stdef', 'blk4_raw_trans1stdef'};
blocksin{15} = {'HollyData_FTBlk1_raw_trans1stdef', 'ft_blk2_raw_trans1stdef', 'ft_blk3_raw_trans1stdef', 'ft_blk4_raw_trans1stdef'};   % as named after Maxfilter
blocksin{19} = {'FTD_blk1_raw_trans1stdef', 'FTD_blk2_raw_trans1stdef', 'FTD_blk3_raw_trans1stdef', 'FTD_blk4_raw_trans1stdef'};   % as named after Maxfilter


%% Define EVENT INFO

eog_thr = [100e-6];     % EOG artefact THRESHOLD (in Volts)

% for stimulus:

epoch  = [-500 2000];    % EPOCH for averaging (milliseconds)
con_labels = {'SPEC_ALL', 'FREE_ALL'};    % Labels of conditions corresponding to trigger codes

offset = zeros(1,2);             % OFFSET between trigger and stimulus presentation, e.g. projector delay (milliseconds)
% for stimulus:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The rest should run smoothly...   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd(bwd) % change to working directory

nr_sbjs = length(subjects); % number of subjects

%Ncons      = length(con_values);
%for c=1:Ncons
%   con_trigs{c} = 'STI101_up'; % look for rising signal in trigger channel for event times
%end
%% open matlabpool if required
%ParType = 0;  % Fun on Login machines (not generally advised!)
%ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
ParType = 2;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)


if ParType
    if matlabpool('size')==0;
        %MaxNsubs = 1;
        %if ParType == 2
        %    for g=1:length(cbu_codes)
        %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
        %    end
        %end
        P = cbupool(19);
        matlabpool(P);
    end
end

load('/imaging/hp02/finger_tapping08/maxFilt/MF-rik/all_bad_channels.mat');
load('/imaging/hp02/finger_tapping08/analysis_spm/preprocessing/new_all_subjs/stimulus/ddm/chan_file');
parfor ss = 1:nr_sbjs
    subjfolder = sprintf('s%d',sname(ss));
    if ~exist( subjfolder, 'dir'); mkdir(subjfolder); end
    cd(subjfolder)
    
    delete *.mat
    delete *.dat
    
    efile = {};
    nr_sess = length( blocksin{ss} );
    
    %clear S D;
    
    swd = fullfile(bwd,subjects{ss});
    fprintf(1, 'Subject: %s\n', swd);
    %try
    %    eval(sprintf('!mkdir %s',swd))
    %end
    %cd(swd)
    
    for ses = 1:nr_sess
        %clear refs
        
        rawfile  = fullfile(bwd,subjects{ss},sprintf('%s.fif',blocksin{ss}{ses}));
        
        fprintf(1, 'Processing %s\n', rawfile);
        
        %trl_name = sprintf('../maxFilt/converted/%s/S_blk%d.mat',sname{ss}, ses);
        %load(trl_name);
        %% Convert:
        S = [];
        S.dataset = rawfile;
        %tmp = load(chanfile);
        S.channels = label;%tmp.label;
        S.outfile  = sprintf('s%d_%s', sname(ss), blkout{ses}) %_%s',filestem,maxflag)
        D = spm_eeg_convert(S);
        
        
        
        % for some reason fiducials disappear in blocks later on... use
        % this to put them back in again:
        if  ses ==1
            blk1_fid = D.fiducials;
        else
            D = fiducials(D, blk1_fid);
        end
        
        %ch = find(strcmp(D.chanlabels,'STI101'));
        %D(ch,:,:) = D(ch,:,:)/10^7;   % Downscale trigger channel so easier to see EOG when display "Other" channels
        
        
        % Correct the trial values in D to correct for inversion of fourth
        % button response:
        %         filename = D.fname;
        %         correct_trial_values(subjfolder, filename, sname(ss), ses);
        %
        %         bin_trials(subjfolder, filename, sname(ss), ses, ss);
        %load(filename)
        
        
        %% Prepare channels 
        if ss <14
            S = [];
            S.D = D;
            S.task = 'settype';
            S.ind = [367 368];
            S.type = 'EOG';
            S.save = 1;
            D = spm_eeg_prep(S);
        end
        
        
        
        %% Correct the trial values in D to correct for inversion of fourth
        % button response:
        filename = D.fname;

        correct_trial_values_new(subjfolder, filename, sname(ss), ses);

        %load(filename)
        disp('corrected');
        
        %% Reference to the average: Diff in 12, need to change
        
        S=[]; S.D = D.fname;
        S.refchan = 'average';
        D = spm_eeg_reref_eeg(S);
        
        
        %
        %% Notch filter
        S = [];
        S.D = D.fname;
        S.type = 'butterworth';
        S.band = 'stop';
        S.freq = [49 51];
        S.dir = 'twopass';
        S.order = 5;
        S.prefix = 'f';
        D = spm_eeg_filter(S);
        
        
        
        %% Notch filter (120 only!)
        S = [];
        S.D = D.fname;
        S.type = 'butterworth';
        S.band = 'stop';
        S.freq = [99 101];
        S.dir = 'twopass';
        S.order = 5;
        S.prefix = 'f';
        D = spm_eeg_filter(S);
        
        
        
        %% Low pass filter
        %         S = [];
        %         S.D = D.fname;
        %         S.type = 'butterworth';
        %         S.band = 'low';
        %         S.freq = 120;
        %         S.dir = 'twopass';
        %         S.order = 5;
        %         S.prefix = 'f';
        %         D = spm_eeg_filter(S);
        %
        
        
        %% high pass filter
        S = [];
        S.D = D.fname;
        S.type = 'butterworth';
        S.band = 'high';
        S.freq = 0.1;
        S.dir = 'twopass';
        S.order = 5;
        S.prefix = 'f';
        D = spm_eeg_filter(S);
        
        
        
        
        %% Epoching:
        disp('epoching')
        %spm('defaults', 'eeg');
        
        S = [];
        S.D = D.fname;
        S.timewin = [-500 2000];
        S.trialdef(1).conditionlabel = 'SPEC_ALL';
        S.trialdef(1).eventvalue = [1:4];
        S.trialdef(2).conditionlabel = 'FREE_ALL';
        S.trialdef(2).eventvalue = [5];
        S.trialdef(3).conditionlabel = 'NULL';
        S.trialdef(3).eventvalue = [6];
     
        
        
        for c=1:2
            S.trialdef(c).eventtype = 'STI101_up'; % look for rising signal in trigger channel for event times
            S.trialdef(c).trlshift = 0;
        end
        
        S.bc = 0;
        S.prefix = 'e';
        S.eventpadding = 0;
        D = spm_eeg_epochs(S);
        
        
        %% Baseline correction
        S = [];
        S.D = D.fname;
        S.timewin = [-500; 0];
        S.prefix = 'b';
        D = spm_eeg_bc(S);
        
        
        %% efile recording
        efile{ses} = D.fname;
        
    end % of ses loop

    
    %% Concatenation of sessions
    %for blk = 1:4
%         efile1 = sprintf('befffMs%d_FT_blk1.mat',sname(ss));
%         efile2 = sprintf('befffMs%d_FT_blk2.mat',sname(ss));
%         efile3 = sprintf('befffMs%d_FT_blk3.mat',sname(ss));
%         efile4 = sprintf('befffMs%d_FT_blk4.mat',sname(ss));
        %efile2{blk} = sprintf('befffMs%d_FT_blk%d.mat',sname(ss), blk);
    %end
    
    S=[];
    S.D = strvcat(efile);%(efile1, efile2, efile3, efile4);
    S.recode = 'same';
    D = spm_eeg_merge(S);
    

    %% Remove Bad channels observed in experiment and maxfilter
    
    
    D = badchannels(D, indchannel(D,all_bad_channels{ss}), 1);
    parsave(D.fname, D);
    
    %% Artifact rejection
    
    S = [];
    S.D = D.fname%sprintf('cbefffMs%d_FT_blk1.mat', sname(ss));%D.fname;
    S.mode = 'reject';
    S.badchanthresh = 0.2;
    
    %
    S.methods(1).channels = {'EOG'};
    S.methods(1).fun = 'peak2peak';
    %S.methods(1).fun = 'threshchan';
    S.methods(1).settings.threshold = 200;
    
    
    S.methods(end).channels = {'MEGPLANAR'};
    S.methods(end+1).fun = 'peak2peak';
    S.methods(end).settings.threshold = 200; %900ft%meg_thr(ss);
    %
    S.methods(end).channels = {'EEG'};
    S.methods(end+1).fun = 'peak2peak';
    S.methods(end).settings.threshold = 200;
    
    D = spm_eeg_artefact(S);
    
    nbadchan(ss) = length(D.badchannels);
    nrej = 0;
    
    disp('done');
    
    %% MEGCOMB
    
    S = [];
    S.D = D.fname;
    S.mode = 'replace';
    D = spm_eeg_combineplanar(S);
    
    
    %     %% Robust averaging
    %     S = [];
    %     S.D = D.fname;
    %     S.robust.ks = 3;
    %     S.robust.bycondition = true;
    %     S.robust.savew = false;
    %     S.robust.removebad = false;
    %     S.circularise = false;
    %     S.prefix = 'm';
    %     D = spm_eeg_average(S);
    %
    %
    %
    %
    %     %% Low pass filter (again for robust averaging)
    %     S = [];
    %     S.D = D.fname;
    %     S.type = 'butterworth';
    %     S.band = 'low';
    %     S.freq = 250;
    %     S.dir = 'twopass';
    %     S.order = 5;
    %     S.prefix = 'f';
    %     D = spm_eeg_filter(S);
    %
    %     %% Holly added: Downsampling:
    %
    %     S = [];
    %     S.D = D.fname;
    %     S.fsample_new = 250;
    %     D = spm_eeg_downsample(S);
    %
    
    %% Delete files to keep space usage down
    %delete s*
    delete Ms*
    delete fM*
    delete ff*
    delete ef*
    delete cb*
    delete mis*
    delete s*
    delete be**
    
    
    
    
    %% Compute Contrasts
    %     S = [];
    %     S.D = D.fname;%sprintf('/imaging/hp02/finger_tapping08/analysis_spm/preprocessing/stimulus/1-40hz/s%d/dfmacbefffMs%d_FT_blk1.mat', sname(ss), sname(ss));
    %     S.c = [1 -1];
    %     S.label = {
    %
    %         'spec-free'
    %
    %         }';
    %     S.WeightAve = 1;
    %     D = spm_eeg_weight_epochs(S);
    %
    %     filename = D.fname;
    %     add_chans_types(filename, sname(ss), blkout{1});
    %
    %
    
    cd ../
end % of subjects loop
save batch_params  nbadchan %nrejects %nevents

if ParType
    matlabpool close
end

return