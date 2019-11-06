%% LBA model trial data
% Use the extracted gradients from the data to make models of the data,
% using fminsearch to find the best fit
addpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/ERPs/variable_ndt');
addpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/MEGCOMB');
clc
%clear all
%load('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/ERPs/average_variable_ndt/MEGCOMB/trialdata.mat');



sname = [ 23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];
ROInum = 96;

%% Parfor
ParType =0;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)


if ParType
    if matlabpool('size')==0;
        %MaxNsubs = 1;
        %if ParType == 2
        %    for g=1:length(cbu_codes)
        %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
        %    end
        %end
        P = cbupool(18);%nr_sbjs);
        matlabpool(P);
        
    end
end
%% LBA modelling
% LBArun = 0;
% if LBArun == 1
% lba_stats = cell(18,102);
% for ss = 1:length(sname)
%     ss
%     for ch = 1:102
%     ch
%        
%         lba_stats{ss,ch} = model_fit4Holly_data_alltrials(trialdata,ss,ch);
%         
%         
%     end
%     
% end
% end
%% STATS

% MEAN NDT

for ss = 1:length(sname)
    ndt_samp = trialdata{ss,1}.ndt_samp;
    for roi = 1:ROInum
        all_ndt_split(ss,roi) = lba_stats{ss,roi}.x2(1,1);
        proportion_ndt_split(ss,roi) = all_ndt_split(ss,roi)/ndt_samp;
    end
end

mean_ndt_split = mean(all_ndt_split,1)*10;
mean_prop_ndt_split = mean(proportion_ndt_split,1);

%%
%% CORRELATIONS

% all time (0-RTs), all trials together
% for ss = 1:length(sname)
%     for roi = 1:ROInum
%         rvalues.allt_alltr_R(ss,roi) = lba_stats{ss,roi}.R;
%         rvalues.allt_alltr_p(ss,roi) = lba_stats{ss,roi}.p;
%         if lba_stats{ss,roi}.p<(0.05/96/17)
%             rvalues.allt_alltr_h(ss,roi) = 1;
%         else
%             rvalues.allt_alltr_h(ss,roi) = 0;
%         end
%     end
% end

% Accumulation period only, all trials together
for roi = 1:ROInum
    roi
    for ss = 1:length(sname)
        ss
        ndt = floor(trialdata{ss,roi}.ndt_samp);
        t1 = round(all_ndt_split(ss, roi));
        if t1 ==0; t1 = 1;end
        if t1>ndt; t1 = ndt; end
        % start arrays to store model and data timeseries across all trials
        clear fullmodel_alltrials fulldata_alltrials
        fullmodel_alltrials = []; % all time
        fulldata_alltrials  = [];
        clear acc_model_alltrials acc_data_alltrials
        acc_model_alltrials = []; % accumulation period only
        acc_data_alltrials  = [];
        cnt = 0;
        for tr = 1:length(trialdata{ss,roi}.trial)
            if ~isempty(trialdata{ss,roi}.trial{1,tr})
                cnt = cnt + 1;
                % keep in samples:
            RT = trialdata{ss,roi}.trial{1,tr}.RT_samp;
            %b0 = trialdata{ss,ch}.trial{1,tr}.b0;
            b1 = trialdata{ss,roi}.trial{1,tr}.LBA_grad;
            t2 = RT - (ndt-t1);
            acc_time = RT-ndt;
            thres = b1*acc_time;
            
            % individual trials:
            
            % all time extraction and correlation
            clear fullmodel fulldata r p
            fullmodel = [zeros(1,t1-1) 0:b1:thres ones(1,RT-t2)*(thres)];
            fulldata = trialdata{ss,roi}.trial{1,tr}.trialdata;
            [r,p] = corrcoef(fullmodel(1:length(fulldata)), fulldata);
            rvalues.allt_r(ss,roi,tr) = r(1,2); rvalues.allt_p(ss,roi,tr) = p(1,2);
            if p(1,2)<(0.05/96/17); rvalues.allt_h(ss,roi,tr) = 1;
            else rvalues.allt_h(ss,roi,tr) = 0; end
            
            % accumulation only extraction and correlation
            clear acc_model acc_data r p 
            acc_model = 0:b1:thres;
            acc_data = trialdata{ss,roi}.trial{1,tr}.trialdata(t1:t2);
            [r,p] = corrcoef(acc_model, acc_data);
            rvalues.acct_r(ss,roi,tr) = r(1,2); rvalues.acct_p(ss,roi,tr) = p(1,2);
            if p(1,2)<(0.05/96/17); rvalues.actt_h(ss,roi,tr) = 1;
            else rvalues.acct_h(ss,roi,tr) = 0; end
            
            % all trials together
            
            fullmodel_alltrials = [fullmodel_alltrials fullmodel]; % all time
            fulldata_alltrials  = [fulldata_alltrials fulldata];
            acc_model_alltrials = [acc_model_alltrials acc_model]; % accumulation period only
            acc_data_alltrials  = [acc_data_alltrials acc_data];
            end
        end
        
        % All trials together correlations:
        clear r p
        [r,p] = corrcoef(fullmodel_alltrials, fulldata_alltrials);
        rvalues.allt_alltrials2_r(ss,roi) = r(1,2); rvalues.allt_alltrials2_p(ss,roi) = p(1,2);
        if p(1,2)<(0.05/96/17); rvalues.allt_alltrials2_h(ss,roi) = 1;
        else rvalues.allt_alltrials2_h(ss,roi) = 0; end
        
        clear r p
        [r,p] = corrcoef(acc_model_alltrials, acc_data_alltrials);
        rvalues.acct_alltrials_r(ss,roi) = r(1,2); rvalues.acct_alltrials_p(ss,roi) = p(1,2);
        if p(1,2)<(0.05/96/17); rvalues.acct_alltrials_h(ss,roi) = 1;
        else rvalues.acct_alltrials_h(ss,roi) = 0; end
    end
end

save('rvalues_pearson.mat', 'rvalues');
