clear all
clc
addpath(genpath('/imaging/hp02/spm12b'));
addpath(genpath('/imaging/local/software/mne'));
addpath(genpath('/imaging/hp02/mmn_08/analysis_spm/new_spm_functions'));
% Root directory for EMEG data
bwd = '/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/';
addpath(genpath(bwd));
%load trial_onsets_indiv_trials.mat

%% open matlabpool if required
%ParType = 0;  % Fun on Login machines (not generally advised!)
%ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
ParType =2;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)


if ParType
    if matlabpool('size')==0;
        %MaxNsubs = 1;
        %if ParType == 2
        %    for g=1:length(cbu_codes)
        %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
        %    end
        %end
        P = cbupool(18);
        matlabpool(P);
    end
end
%% Define loops

% % the inversion methods to be performed
inver_meth_path = { ''};%, 'MEGgrad_MSP/', 'MEGgrad_sLoreta/', 'MEGgrad_Beamf/' };%, 'MEGgrad_Beam/'};
inv_meth        = { 'IID'};%,               'GS',          'LOR', 'EBB' };
time_wind_path  = { 'all_time'};
windows         = { [-500 1500]};

% does preprocessed data need to be copied? 1 = yes, 0 = no:
preproc_copyflag = 0;
sourceloc_path = '/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/';

%% Define subjects' preprocessed files

prefix = 'devracbefffM';
sname = [19 23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];

for i = 1:length(sname)
    if i ==14
    subjects{i} = sprintf('%ss%d_FT_blk1.mat', prefix, sname(i));   % letter strings for subject-specific paths
    subj_dotdat{i} = sprintf('%ss%d_FT_blk1.dat', prefix, sname(i));
    else
        subjects{i} = sprintf('%ss%d_FT_blk1.mat', prefix, sname(i));   % letter strings for subject-specific paths
    subj_dotdat{i} = sprintf('%ss%d_FT_blk1.dat', prefix, sname(i));
    
    end
    MRI_dir{i} = sprintf('/imaging/hp02/mmn_08/mri_images/spm12_mri/s%d/', sname(i));  % directory with this participant's MRI (*.img)
    
end

nr_subs = length(sname);
fprintf(1, 'Going to process %d data sets\n', nr_subs);


%% Define Analysis Parameters/Processing options/ Inversion options that don't change:
% Names of conditions after averaging
% inv_trls = {'SPEC_ALL', 'FREE_ALL', 'NULL', 'ACTION'};    % Labels of conditions corresponding to trigger codes    % Labels of conditions corresponding to trigger codes

freq_start = 0.1;
freq_end = 100;

% Forward model options, for_typ{EEG/MEG}
for_typ{1} = 'EEG BEM'; % EEG
for_typ{2} = 'Single Shell'; % MEG (mags+grads)

% Processing options
segment_only = 0;    % Flag whether only MRI segmentation (up to, not including, coregistration) shall be done (1: only segment; 0: do it all)
redoreg_flag = 0;    % Flag whether to redo MRI and registration model (eg if want to change fids)
redoformod_flag = 1; % Flag whether to redo forward modelling (eg if only one for_typ above)
display_flag = 1;    % display results on the fly
write3Dflag = 1;     % write SPM volumes of results

%Inversion Options

% Which sensor types to use
inv_mods = {'MEGPLANAR'};   % Just MEG gradiometers  %inv_mods = {'MEG';'MEGPLANAR';'EEG'};   % Fusion of mulitple modalities

% Cortical surface option
mesh_size = 2;   % [1-3] for coarse-normal-fine

% Whether to use headshape for coregistration (in datareg)
use_headshape = 1;



%% Loop around each inversion method and each time window:

%for inv_cnt = 1:length(inv_meth)
inv_cnt = 1;


% Define Analysis Parameters
% Which inversion method to use
inv_typ = inv_meth{inv_cnt};   % Minimum Norm Least-Squares %inv_typ = 'GS';   % (MSP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Source localisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear D
forward_modeling = 0;
if forward_modeling == 1
parfor ss = 2:nr_subs
    
    
    % Processing options
    segment_only = 0;    % Flag whether only MRI segmentation (up to, not including, coregistration) shall be done (1: only segment; 0: do it all)
    redoreg_flag = 0;    % Flag whether to redo MRI and registration model (eg if want to change fids)
    redoformod_flag = 1; % Flag whether to redo forward modelling (eg if only one for_typ above)
    display_flag = 1;    % display results on the fly
    
    
    % Make sure in right starting place:
    cd(sourceloc_path)
    
    if ~exist(inver_meth_path{inv_cnt}, 'dir'); mkdir(inver_meth_path{inv_cnt}); end
    %cd(inver_meth_path{inv_cnt}) % inversion method folder
    %
    subjfolder = sprintf('s%d',sname(ss));
    if ~exist( subjfolder, 'dir'); mkdir(subjfolder); end
    cd(subjfolder)
    
    
    
    % current subject's data directory
    subj_dir = pwd;
    fprintf(1, '\n Current subject directory: %s\n\n', subj_dir);
    
    data_path = pwd;%[path_file, inver_meth_path{inv_cnt}, time_wind_path{wind_cnt}, data_name{ss}];
    data_subj = subjects{ss};
    data_file = [data_path,'/', data_subj];
    
    D = spm_eeg_load( data_file );
    
    % Initialise... (if want to be safe!)
    D.inv = {struct('mesh', [], 'gainmat', [])};
    
    current_formod = 1;
    
    indF = [];
    val = 0;
    
    %% MRI processing and coregistration
    if redoreg_flag == 1 || ~isfield(D.inv{1},'mesh') || ~isfield(D.inv{1},'datareg')
        
        val = val + 1;
        D.val = val;
        D.inv{val}.date    = strvcat(date,datestr(now,15));
        D.inv{val}.comment = {sprintf('%s_%s',inv_typ,char(inv_mods)')};    % remember inversion type and sensor configuration
        D.inv{val}.mesh    = [];
        D.inv{val}.datareg = [];
        D.inv{val}.forward = [];
        
        % Locations/names of MRIs
        if sname(ss) ~= 24 %If not subject 24 who doesn't have an MRI MPRAGE
            % Check whether MRI present (care if more than one such *img!)
            D.inv{val}.mesh.sMRI = spm_select('FPList',MRI_dir{ss},'^MPRAGE.*\.img$');
            if isempty(D.inv{val}.mesh.sMRI)
                error('MRI not found')
            end
            
            % s24:
        else
            % Check whether MRI present (care if more than one such *img!)
            D.inv{val}.mesh.sMRI = spm_select('FPList',MRI_dir{ss},'^single_subj_T1.*\.nii$');
            if isempty(D.inv{val}.mesh.sMRI)
                error('MRI not found')
            end
        end
        
        
        %% Check whether inverse-normalised surfaces present
        D.inv{val}.mesh.def  = spm_select('FPList',MRI_dir{ss},'^y_.*\.nii$');
        if isempty(D.inv{val}.mesh.def)
            warning('No inverse normalisation parameters found - segmenting may take a while!')
        end
        
        %% Normalise sMRI (if not done already), and create inverse-normalised surfaces
        D.inv{val}.mesh = spm_eeg_inv_mesh(D.inv{val}.mesh.sMRI, mesh_size);
        if display_flag, spm_eeg_inv_checkmeshes(D); end
        % D.save;
        
        % If coregistration not yet done...
        if segment_only,
            continue;   % Stop processing for this subject here, continue with next one
        end;
        
        % If coregistration done manually and saved...
        newmrifid           = [];
        newmrifid.fid.pnt   = D.inv{val}.mesh.fid.fid.pnt(1:3,:);  % I think these are the fids after "Save"...
        newmrifid.fid.label = {'Nasion';'LPA';'RPA'};
        newmrifid.pnt       = D.inv{val}.mesh.fid.pnt;        % Scalp mesh points from MRI above
        
        meegfid = D.fiducials;
        %% Remove nose points (assume those for which y>0 and z<0)
        meegfid.pnt(find(meegfid.pnt(:,2)>0 & meegfid.pnt(:,3)<0),:) = [];
        
        fprintf(1, 'Coregistering\n');
        D = spm_eeg_inv_datareg_noui_130614(D, val, meegfid, newmrifid,use_headshape);
        
        fprintf(1, 'Done Coregistering\n');
        
        for ind = 1:length(D.inv{val}.datareg)
            fprintf(1, '%d ', ind);
            d = D.inv{val}.datareg(ind).fid_eeg.fid.pnt - D.inv{val}.datareg(ind).fid_mri.fid.pnt;
            %err(ss,ind,val) = mean(sqrt(sum(d.^2,2)));
        end
        fprintf(1, '\n');
        
        redoreg_flag = 0;
        
    else
        D.inv{val}.mesh    = D.inv{1}.mesh;
        D.inv{val}.datareg = D.inv{1}.datareg;
    end % if redoreg_flag...
    
    %% Computing forward model/leadfield
    if redoformod_flag == 1 | ~isfield(D.inv{1},'forward') | ~exist(D.inv{1}.gainmat)
        
        fprintf(1, 'Creating forward model\n');
        %% Create forward model (BEM) (could conditionalise this bit on modality inverted...)
        
        D.inv{val}.forward = struct([]);
        for ind = 1:length(D.inv{val}.datareg)
            D.inv{val}.forward(ind).voltype = for_typ{ind};
        end
        
        fprintf(1, 'Computing leadfield\n');
        D = spm_eeg_inv_forward(D);
        
        if display_flag
            for ind = 1:length(D.inv{val}.datareg)
                spm_eeg_inv_checkforward_hp_130614(D, val, ind); %pause(3);
            end
        end
        current_formod = val;
        redoformod_flag = 0;
        D.save;
    else
        D.inv{val}.forward = D.inv{current_formod}.forward;
        D.inv{val}.gainmat = D.inv{current_formod}.gainmat;
        D.save;
    end     % redoformod_flag = ...
end  
end

clear D
%get the trial names ready for inversion:
% clear inv_trls
% for ss = 1:nr_subs
%     for i  = 1:size(tr_onsets{ss},1)
%         inv_trls{ss}{i} = tr_onsets{ss}{i,1};
%     end
% end    
%% Invert & Contrast
inv_trls = {};
for ss = 2:nr_subs
    data_path = sprintf('%s/s%d',bwd, sname(ss));%[path_file, inver_meth_path{inv_cnt}, time_wind_path{wind_cnt}, data_name{ss}];
    data_subj = subjects{ss};
    data_file = [data_path,'/', data_subj];
    
% extract trial names
    load(data_file);
    data_fname = D.data.fname;
    
    for tr = 1:length(D.trials)
        if ~strcmp(D.trials(1,tr).label, 'NULL')
            inv_trls{ss}{tr} = D.trials(1,tr).label;
        end
    end
end
    
for ss = 2:nr_subs    
    ss
    write3Dflag = 0;     % write SPM volumes of results
    
    
    
    % Load forward modelled file:
    cd(sourceloc_path)
    %cd(inver_meth_path{inv_cnt}) % inversion method folder
    subjfolder = sprintf('s%d',sname(ss));
    cd(subjfolder)

    % current subject's data directory
    data_path = pwd;%[path_file, inver_meth_path{inv_cnt}, time_wind_path{wind_cnt}, data_name{ss}];
    data_subj = subjects{ss};
    data_file = [data_path,'/', data_subj];
    
    
    D = spm_eeg_load( data_file );
    %D.fnamedat = data_fname;
    
    val = length(D.inv);
    
    D.inv{val}.inverse = [];   % Clear to be safe!
    D.inv{val}.inverse.trials = inv_trls{ss};
    D.inv{val}.inverse.type   = inv_typ;
    
    D.inv{val}.inverse.woi    = [-500 1500];
    D.inv{val}.inverse.lpf    = freq_start;
    D.inv{val}.inverse.hpf    = freq_end;
    
    D.inv{val}.inverse.modality = inv_mods;
    
    D = spm_eeg_invert(D);
    D.save;
    for wind_cnt = 1:length(time_wind_path)
        D.inv{val}.contrast.woi  = [windows{wind_cnt}(1) windows{wind_cnt}(2)];
        %D.inv{val}.contrast.fboi = [freq_start freq_end];
        D.inv{val}.contrast.type = 'evoked';
        D2 = D;
        D2 = spm_eeg_inv_results(D2);
        
        %% Create images for stats with smoothing:
        
        
        % Write result to SPM volumes
        if write3Dflag
            D2.inv{val}.contrast.smoothing = 8;
            
            D2 = spm_eeg_inv_Mesh2Voxels(D2);
            %SourceImgs{val}{ss} = strvcat(D.inv{val}.contrast.fname);
        end
        
        %indF(ss,val) = D.inv{val}.inverse.F;
        
        
    end
    
    
    
    
end

%end
if ParType
    matlabpool close force CBU_Cluster
end