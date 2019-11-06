clc
clear all

cwd = '/imaging/hp02/finger_tapping08/analysis_spm/time_freq_analysis/new_subs/stimulus/gamma/rm_rtf_old';
%spm('defaults', 'eeg');

sname = [18 19 23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];
prefix1 = 'rmbad_rtf_acbefffM';

subjects= cell(length(sname), 1);
subj_dotdat = cell(length(sname), 1);

for ss = 1:length(sname)
    %data_name{ss} = sprintf('/s%d/', sname(ss));
    subjects{ss} = sprintf('%s/%ss%d_FT_blk1.mat', cwd,  prefix1, sname(ss));
    subj_dotdat{ss} = sprintf('%s/%ss%d_FT_blk1.dat', cwd,  prefix1, sname(ss));
    
    
    %     S = [];
    %     S.D = subjects{ss};
    %     S.prefix = 'rmbad_';
    %     D = spm_eeg_remove_bad_trials(S);
end

%%

clear good_RT
disp('Calculate RTs')

good_RT = cell(1,length(sname));

for ss = 1:length(sname)
    clear D
    ss;
    load(subjects{ss});
    for i = 1:length(D.trials)
        j = 1;
        stim = 0;
        resp = 0;
        while j<=length(D.trials(1,i).events(:,1)) && resp ==0
            if strcmp(D.trials(1,i).events(j,1).type, 'STI101_up') && stim == 0
                stim_time = D.trials(1,i).events(j,1).time;
                stim = 1;
            elseif strcmp(D.trials(1,i).events(j,1).type, 'STI101_up') && stim == 1
                resp_time = D.trials(1,i).events(j,1).time;
                resp = 1;
                good_RT{ss}.RTs(i,1) = resp_time-stim_time;
                good_RT{ss}.cond{i,1} = D.trials(1,i).label;
                if strcmp(good_RT{ss}.cond{i,1}, 'SPEC_ALL')
                    good_RT{ss}.condN{i,1} = 1;
                elseif strcmp(good_RT{ss}.cond{i,1}, 'FREE_ALL')
                    good_RT{ss}.condN{i,1} = 2;
                end
                good_RT{ss}.finger(i,1) = D.trials(1,i).events(j,1).value;
            end
            
            j= j+1;
        end
    end
end
%% Calculate Single trial Drifts

% Output from LBA, copied over from: /imaging/hp02/finger_tapping08/behaviour/RTs/scripts_ch4_HOLLY/scripts/all_free_params
load('/imaging/hp02/finger_tapping08/behaviour/RTs/scripts_ch4_HOLLY/scripts/all_free_params/par_Free_22.mat')
load('/imaging/hp02/finger_tapping08/behaviour/RTs/scripts_ch4_HOLLY/scripts/all_free_params/par_Spec_22.mat')

% Initialise cells
EAA_tot = cell(length(sname), 1);
EAG_tot = cell(length(sname), 1);
report = cell(length(sname), 2);

for ss = 1:length(sname)
    sas = ss;
    % To account for the 2 ECoG patients in the middle of this dataset:
    if ss>14
        sas = ss+2;
    elseif ss > 17 % to take into account the removable of subject s532
        sas = ss+3;
        
    end
    
    b = par_Free(sas,1).B(1); % Same across fingers and conditions
    vm(1,1:4) = par_Spec(sas,1).Ame; % Mean drift rate for spec - different for different conditions and fingers!
    vm(2,1:4) = par_Free(sas,1).Ame;
    vs(1,1:4) = par_Spec(sas,1).Astd; % Mean drift rate for spec - different for different fingers!
    vs(2,1:4) = par_Free(sas,1).Astd;
    C0d2 = (par_Free(sas,1).C0(1))/2; % max bias, fixed across all
    t0 = par_Free(sas,1).T0; % NDT, same across conds
    
    cnt = 0;
    % each trial
    for i = 1:length(good_RT{1,ss}.RTs(:,1))
        clear EAAw EAAL EAGw EAGL ti finger cond
        ti = good_RT{1,ss}.RTs(i,1);
        
        cond = good_RT{ss}.condN{i,1}; % Spec = 1, Free = 2
        
        if cond ==1
            finger = good_RT{ss}.finger(i,1)-10; % Which finger was pressed
        elseif cond ==2
            finger = good_RT{ss}.finger(i,1)-50;
        end
        
        % EEA of the winning accumulator
        EAAw = 0.5*(b-C0d2)*(ti-t0); % need to check whether + or -
        Ev_w = (b-C0d2)/(ti-t0);
        EAGw = (b-C0d2)/(ti-t0);
        
        % other fingers (losers)
        all_fingers = 1:4;
        losing_fingers = all_fingers(all_fingers~=finger); % just get the losing fingers
        Ev_l = zeros(1,3);
        EAAL= zeros(1,3);
        EAGL= zeros(1,3);
        for f = 1:3
            if cond == 2
                vmL = vm(2,losing_fingers(f)); % mean acc of the loser
                vsL = vs(2,losing_fingers(f)); % std acc of the loser
                
                x = (Ev_w-vmL)/vsL;
                phi = (1/sqrt(2*pi)) * exp(-(x^2)/2);
                
                fun = @(xf) exp(-(xf.^2)./2);
                
                PHI = (1/sqrt(2*pi)) * integral(fun, -Inf, x);
                
                Ev_l(f) = vmL - vsL*(phi/PHI);
                
                EAAL(f) = 0.5*((Ev_l(f)*(ti-t0))) * (ti-t0); % I am unsure about the +
                EAGL(f) = (Ev_l(f)*(ti-t0))/(ti-t0);
                
                % set NaNs to zero
                if isnan(EAAL(f))
                    EAAL(f) = 0;
                end
                if isnan(EAGL(f))
                    EAGL(f) = 0;
                end
                
                
                
            end
            
            
            
        end
        
        % if the reaction time is less than the non decision time,
        % something is odd in their response and this trial should
        % be removed from further analysis. Set EAA and EAG to
        % zero and record this trial. I will have to remove these
        % trials manually later, before the correlation.
        if ti < t0
            EAAL = zeros(1,3);
            EAGL = zeros(1,3);
            EAAw = 0;
            EAGw = 0;
            cnt  = cnt+ 1;
            report{ss,1}(cnt) = i;
            report{ss,2} = cnt;
        end
        
        % Sum up all the EAAs and EAGs (winner and losers)
        if cond == 1
            EAA_tot{ss}(i) = EAAw;
            EAG_tot{ss}(i) = EAGw;
        elseif cond == 2
            EAA_tot{ss}(i) = EAAw + sum(EAAL);
            EAG_tot{ss}(i) = EAGw + sum(EAGL);
        end
        
        
        
        
    end
    
    
end

for ss = 1:length(sname)
    %drift_trials(ss,1) = length(drift_rate{ss,1}.rates);
    
end

%% Save drifts

save('ndt_violation_report.mat', 'report');
save('EAA.mat', 'EAA_tot');
save('EAG.mat', 'EAG_tot');








