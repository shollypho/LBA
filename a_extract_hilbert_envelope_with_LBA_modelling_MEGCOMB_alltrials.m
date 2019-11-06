%clear all
clc

addpath(genpath('/imaging/hp02/spm12b'));
addpath(genpath('/imaging/local/mne'));
addpath(genpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/ERPs/variable_ndt'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/parcellations_AleT_100Hzds_m500_1500ms');
%% Subjects, LBA parameters and channel names
sname = [23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];
%path_rtf_files = '/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/parcellations';
prefix_temp = 'tempROI_devracbefffMs';
prefix = 'devracbefffMs';
%s_idx = [2:12, 15:18, 20,21]; % corresponding LBA parameters as some subjects are excluded

% Load the LBA parameters
load('/imaging/hp02/finger_tapping08/behaviour/RTs/scripts_ch4_HOLLY/scripts/all_free_params/par_Free_22.mat')
load('/imaging/hp02/finger_tapping08/behaviour/RTs/scripts_ch4_HOLLY/scripts/all_free_params/par_Spec_22.mat')

% Load mask data
%load mask_1stEigen.mat

HO_list;

ROInum=96;

%% open matlabpool if required
% %     ParType = 0;  % Fun on Login machines (not generally advised!)
% %     ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
ParType = 0;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)


if ParType
    if matlabpool('size')==0;
        %MaxNsubs = 1;
        %if ParType == 2
        %    for g=1:length(cbu_codes)
        %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
        %    end
        %end
        P = cbupool(17);
        P.ResourceTemplate='-l nodes=^N^,mem=16GB,walltime=72:00:00';
        matlabpool(P);
    end
end
%% Extract the hilbert Envelope
% frequency windows of interest
fHp = [12 18 30 48 60];
fLp = [30 20 60 60 90];
for band=5
%         band
         fHpi = fHp(band);
         fLpi = fLp(band);
        %
        disp('extract hilbert transform')
        for ss = 1:length(sname)
    
            ss
            data = [];
            for roi = 1:ROInum
                
                for tr = 1:size(mask_1stEigen(ss).data,3)
                    for t = 1:size(mask_1stEigen(ss).data,2)
                        data(roi,t,tr) = mask_1stEigen(ss).data(roi,t,tr);
                    end
                end
            end
            
            fChan{ss,1}  = envelope_ROI(data, fHpi, fLpi);
        end
    %
    
        disp('done')
        cd(sprintf('%d-%dHz', fHpi, fLpi));
        %     save  ('fChanRS.mat','fChanRS','fChanRS_tr', 'fChanmean');
        %
        %% LBA modelling
        % Looping through subjects
        disp('Extracting single trial LBA estimated gradients and data');
        for ss = 1:length(sname)
    
            ss
            
            % get trial number from inversion file
            fname = sprintf('../s%d/%s%d_FT_blk1.mat', sname(ss),prefix_temp,sname(ss));
            
            load(fname)
            tr_num = length(D.trials);
            
            % load original D struct for the events
            fname = sprintf('../s%d/%s%d_FT_blk1.mat', sname(ss),prefix,sname(ss));
            
            load(fname)
            
                
            % Getting the correct index for each subject (a bit of a hack)
            sas = ss+2; % because not including s18 or s19 atm, need to re-maxfilter and preprocess s18.
            if ss>12; sas = ss+4; % To account for the 2 ECoG patients in the middle of this dataset:
            elseif ss > 16;sas = ss+5; % to take into account the removable of subject s532
            end
    
            % Load the LBA parameters for this subject:
            clear bm vm vs C0d2 t0
            b = par_Free(sas,1).B(1); % Same across fingers and conditions
            vm(1,1:4) = par_Spec(sas,1).Ame; % Mean drift rate for spec - different for different conditions and fingers!
            vm(2,1:4) = par_Free(sas,1).Ame; % Free Choice
            vs(1,1:4) = par_Spec(sas,1).Astd; % Mean drift rate for spec - different for different fingers!
            vs(2,1:4) = par_Free(sas,1).Astd; % Free choice
            C0d2 = (par_Free(sas,1).C0(1))/2; % max bias, fixed across all
            t0 = par_Free(sas,1).T0; % NDT, same across conds
    
   
            % counters:
            i = 0; % trials
            i_f = 0; % free trials
            i_s = 0; % specified trials
            cnt = 0;
    
    
            %% Trials
            for tr = 1:tr_num
    
                j = 1; stim = 0; resp = 0;
    
                if strcmp(D.trials(1,tr).label,'NULL')
                    % Null trial, skip
                    test=1;
                else
                    %% Get the REACTION TIME from this trial to calculate the drift rate
                    % loop through the D.trials.events to calculate the RTs
    
                    while j<=length(D.trials(1,tr).events(:,1)) && resp ==0
    
                        % if the first STI101_up, should be the stimulus
                        % presentation, get time
                        if strcmp(D.trials(1,tr).events(j,1).type, 'STI101_up') && stim == 0
    
                            stim_time = D.trials(1,tr).events(j,1).time;
                            stim = 1;
    
                            % if the second STI101_up, should be the button
                            % response, get time and trial type condiction (spec or free)
                        elseif strcmp(D.trials(1,tr).events(j,1).type, 'STI101_up') && stim == 1
    
                            resp_time = D.trials(1,tr).events(j,1).time;
                            resp = 1; % use to break out of while loop
                            good_RT.RTs = resp_time-stim_time; % RT
                            good_RT.cond = D.trials(1,tr).label; % Trial condition
    
                            % assign code for each trial condition
                            if strncmpi(good_RT.cond, 'SPEC_ALL',8)
                                good_RT.condN = 1;
                            elseif strncmpi(good_RT.cond, 'FREE_ALL',8)
                                good_RT.condN = 2;
                            end
    
                            good_RT.finger = D.trials(1,tr).events(j,1).value; % which finger was pressed
                        end
    
                        j= j+1;
                    end
    
                    clear ti
    
                    ti = good_RT.RTs;
    
                    % Make sure to only work with trials that have RTs
                    % greater than the NDT (t0)
                    if ti > t0
                        i = i + 1;
    
                        %% Drift rate - calculated from the LBA parameters and single trial RT
    
                        clear EAAw EAAL EAGw EAGL finger cond
    
    
                        cond = good_RT.condN; % Spec = 1, Free = 2
    
                        if cond ==1
                            finger = good_RT.finger-10; % Which finger was pressed as spec fingers: 11-Index, 12 = Middle, 13 = Ring, 14 = Little
                        elseif cond ==2
                            finger = good_RT.finger-50; % Which finger was pressed as free fingers: 51-Index, 52 = Middle, 53 = Ring, 54 = Little
                        end
    
                        % EEA of the winning accumulator
                        %EAAw = 0.5*(b-C0d2)*(ti-t0); % Winning Area: 1/2 base x height
                        Ev_w = (b-C0d2)/(ti-t0);
                        EAGw = (b-C0d2)/(ti-t0); % Winning gradient: height/base (y/x)
    
                        % Sum up all the EAAs and EAGs (winner and losers)
                        if cond == 1
                            i_s = i_s + 1;
    
                            EAG_tot_spec{ss}(i_s) = EAGw;%EAA_tot{ss}(i) = EAAw;
                            EAG_tot_all{ss}(i) = EAGw;
    
                            % If free choice condition, also need to calculate accumution of the losers:
                        elseif cond == 2
                            i_f = i_f + 1;
    
                            all_fingers = 1:4;
                            losing_fingers = all_fingers(all_fingers~=finger); % just get the losing fingers
                            Ev_l = zeros(1,3);
                            %EAAL= zeros(1,3);
    
                            for f = 1:3 % loop through 3 losing fingers
                                if cond == 2 % if the free choice trial
                                    vmL = vm(2,losing_fingers(f)); % mean acc of the loser
                                    vsL = vs(2,losing_fingers(f)); % std acc of the loser
    
                                    x = (Ev_w-vmL)/vsL;
                                    phi = (1/sqrt(2*pi)) * exp(-(x^2)/2);
    
                                    fun = @(xf) exp(-(xf.^2)./2);
    
                                    PHI = (1/sqrt(2*pi)) * integral(fun, -Inf, x);
    
                                    Ev_l(f) = vmL - vsL*(phi/PHI);
    
                                    EAAL(f) = 0.5*((Ev_l(f)*(ti-t0))) * (ti-t0); % I am unsure about the + c0/2, I think the term should not be there
    
    
                                    % set NaNs to zero, they do not contribute to
                                    % the overall accumulation
                                    if isnan(EAAL(f))
                                        EAAL(f) = 0;
                                        Ev_l(f) = 0;
                                    end
    
    
                                end
    
                            end
    
                            %EAA_tot{ss}(i) = EAAw + sum(EAAL);
                            %new_height = (2*EAA_tot{ss}(i))/(ti-t0); % to calculate the new gradient, need the height of the triangle = 2 x total area/base
    
                            % Sum the gradients over winner and losers
                            EAG_tot_free{ss}(i_f) = EAGw + sum(Ev_l);  %new_height/(ti-t0); % gradient = y/x = h/b = new height/ RT-ndt
                            EAG_tot_all{ss}(i) = EAGw + sum(Ev_l);
                        end
    
                        %% Caluculate the Gradient of the accumulated activity
    
                        samp_time = 1/100; % time of each sample in the data
                        baseline_samp = 0.5/samp_time;
                        RT_samp = round(ti/samp_time);
    
                        %% Channels
                        % Go through each channel, to calculate single trial drift rates and
                        % extract gradients from the channel recordings
                        for roi = 1:ROInum
    
    
                            % load the ERP data - just from 0ms-RTms - so remove baseline and data after RT:
                            trialdata{ss,roi}.trial{tr}.trialdata(1:RT_samp) = fChan{ss}(roi, baseline_samp:RT_samp+baseline_samp-1, tr);
                            trialdata{ss,roi}.trial{tr}.RT = ti;
                            trialdata{ss,roi}.trial{tr}.RT_samp = RT_samp;
                            trialdata{ss,roi}.ndt = t0;
                            trialdata{ss,roi}.ndt_samp = t0/samp_time;
                            trialdata{ss,roi}.trial{tr}.LBA_grad = EAG_tot_all{ss}(i);
                            if cond == 1
                                trialdata{ss,roi}.trial{tr}.LBA_spec = EAG_tot_spec{ss}(i_s);
                                trialdata{ss,roi}.trial{tr}.cond = 1;
                                trialdata{ss,roi}.trial{tr}.triallabel = 'SPEC';
    
                            elseif cond == 2
                                trialdata{ss,roi}.trial{tr}.LBA_free = EAG_tot_free{ss}(i_f);
                                trialdata{ss,roi}.trial{tr}.cond = 1;
                                trialdata{ss,roi}.trial{tr}.triallabel = 'FREE';
                            end
    
                            trialdata{ss,roi}.trial{tr}.b0 = mean(fChan{ss}(roi, baseline_samp-(0.1/samp_time):baseline_samp-1, tr)); % Mean of the 100ms before 0

                        end
    
    
    
                    end % if ti>t0
                end % if not a null trial
            end %trials
    
    
    
        end % ss
    
       save('trialdata.mat', 'trialdata');
    %% Modelling
    disp('Model fit');
    
    
    
    lba_stats = cell(length(sname),ROInum);
    for ss = 1:length(sname)
        ss
        for roi = 1:ROInum
            roi
            
            %for tr = 1:length(trialdata{ss,1}.trial)
            %   if (~isempty(trialdata{ss,ch}.trial{1,tr}))
            
                lba_stats{ss,roi} = corr_model_fit4Holly_data_alltrials(trialdata,ss,roi);
            %  end
           
            %end
            
            
        end
        
    end
    
    
    
    
    
    % save all workspace
    %clear fChan 
    date1 = date;
    %save(sprintf('workspace_%s_fband_%d_%dHz_MegComb.mat', date1, fHpi, fLpi));
    save('lba_stats.mat', 'lba_stats');
    
    cd ../
end % Bands
% if ParType
%     matlabpool close force CBU_Cluster
% end
%% provided lba_stats saves properly, can then do stats.

