function lba_stats=corr_model_fit4Holly_data_alltrials(trialdata, ss, roi)



ndt =round((trialdata{ss,roi}.ndt)/0.01);
%data = trialdata{ss,ch}.trial{tr}.trialdata;
data =[];
cnt = 0;
b0_all = [];b1_all = [];RT_all = [];
for tr = 1:length(trialdata{ss,roi}.trial)
    if (~isempty(trialdata{ss,roi}.trial{1,tr}))
        cnt = cnt+1;
        % extract all the LBA parameters modelled from the behavioural data
        % and the reaction times:
        RT_all(cnt) = trialdata{ss,roi}.trial{tr}.RT_samp;
        b1_all(cnt) = trialdata{ss,roi}.trial{tr}.LBA_grad;
        b0_all(cnt) = double(trialdata{ss,roi}.trial{tr}.b0);
        
        
        % ROI data, concatinated over trials
        data =[data, trialdata{ss,roi}.trial{1,tr}.trialdata];
        
    end
end
delta = @(x,t) x>t;
%x = 1:length(data);
%% Iterate through lags for NDT
cnt = 0;
for lag = 1:1:ndt
    %lag
    regrm_tot = [];
    for tr = 1:length(b0_all)
        lbasig = [];
        b1=  b1_all(tr);
        RT = RT_all(tr);
        
        %t1 pre-t0 t2 post b1 slope/drift
        yt = @(x1,b1,t1,t2) b1.*(delta(x1,t1).*(x1-t1) - delta(x1,t2).*(x1-t2));
        
        t1=lag;%lag(i)
        t2=ndt-t1;
        t2start = RT-(ndt-t1);
        x1 = 1:RT;
        t2 = RT - (ndt-t1);
        
        acc_time = RT-ndt;
        thres = b1*acc_time;
        
        %lbasig = [zeros(1,t1) b1:b1:thres ones(1,RT-t2start)*thres ];
        
        
        lbasig = yt(x1,b1,t1,t2); % using old fmin version
        %ramp = lbasig(t1:end); 
        %regrm_tr = [zeros(1,t1), ramp,ones(1,t2).*ramp(end)];%regrm(i)
        
        regrm_tot = [regrm_tot, lbasig];
    end
    cnt = cnt + 1;
    [rho(cnt),pval(cnt)] = corr(regrm_tot', data','type','Spearman');%rho(i) for each lag
    [rr,pp] = corrcoef(regrm_tot', data');
    r(cnt) = rr(1,2);
end


 
%% LBA Stats

[maxR, t1_opt_idx] = max(abs(rho));

lba_stats.R = maxR;
lba_stats.p = pval(t1_opt_idx);
lba_stats.rho = rho;

t1=t1_opt_idx;%lag(i)
t2=ndt-t1;

lba_stats.x2 = t1;
lba_stats.t1 = t1;
% 
% opt_t1_model = [];
% for tr = 1:length(b0_all)
%     
%     b1=  b1_all(tr);
%     RT = RT_all(tr);
%     
%     x1 = 1:RT;
%     lbasig = yt(x1,b1,t1,t2);
%     
%     opt_t1_model = [opt_t1_model, lbasig];
% end
% 
% %figure; plot(opt_t1_model/20000); hold on; plot(data, 'r')
% 
% [rho_opt,pval_opt] = corr(opt_t1_model', data','type','Spearman');

