%% Stats to compare action trials to Scrambled trials and spec to free trials
clear spec_trl_idx free_trl_idx act_trl_idx

clc

%% Extract indices for specified and for free trials
sname = [23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];
ROInum = 96;

for ss = 1:17
    

    
    act_trl_num(ss) = length(trialdata{ss,1}.trial);
    spec_cnt = 0;
    free_cnt = 0;
    empty_cnt = 0;
    act_cnt = 0;
    for tr = 1:act_trl_num(ss)
        if isfield(trialdata{ss,1}.trial{1,tr}, 'triallabel')
            if strcmp(trialdata{ss,1}.trial{1,tr}.triallabel, 'SPEC')
                spec_cnt = spec_cnt + 1; act_cnt = act_cnt + 1;
                spec_trl_idx{ss}(spec_cnt) = tr;
                act_trl_idx{ss}(act_cnt) = tr;
                
            elseif strcmp(trialdata{ss,1}.trial{1,tr}.triallabel, 'FREE')
                free_cnt = free_cnt + 1;act_cnt = act_cnt + 1;
                free_trl_idx{ss}(free_cnt) = tr;
                act_trl_idx{ss}(act_cnt) = tr;
                
            end
        else
            empty_cnt = empty_cnt + 1;
            empty_tr{ss}(empty_cnt) = tr;
        end
        
    end
end

%% Convert r values to Fisher z values

for ss = 1:length(sname)
    ss
    for ch = 1:ROInum
        s_cnt = 0; f_cnt = 0; % counters for free and spec trials
        %for tr = 1:length(act_trl_idx{ss})
        % all
        r = rvalues.allt_alltrials_r(ss,ch);
        if r == -1; r = -0.9999; elseif r == 1; r = 0.9999;end % fix for r = -1 or r = 1;
        zvalues.allt_alltrials_z{ss,ch} = .5.*log((1+r)./(1-r));
        
        % Scrambled trials (already averaged across trials)
        for permi = 1:10000
            r = rvaluesScram.all_r(permi,ss,ch);
            if r == -1; r = -0.9999; elseif r == 1; r = 0.9999;end % fix for r = -1 or r = 1;
            zvaluesScram.act_all_z(permi,ss,ch) = .5.*log((1+r)./(1-r));
        end
        
        % Spec
        r = rvalues.allt_alltrials_r_spec(ss,ch);
        if r == -1; r = -0.9999; elseif r == 1; r = 0.9999;end % fix for r = -1 or r = 1;
        zvalues.allt_alltrials_z_spec{ss,ch} = .5.*log((1+r)./(1-r));
        
        % Free
        r = rvalues.allt_alltrials_r_free(ss,ch);
        if r == -1; r = -0.9999; elseif r == 1; r = 0.9999;end % fix for r = -1 or r = 1;
        zvalues.allt_alltrials_z_free{ss,ch} = .5.*log((1+r)./(1-r));
        
    end
end

%% Test permutation distributions are normal
for ss = 1:17
    for ch = 1:96
        h(ss,ch) = kstest(zvaluesScram.act_all_z(:,ss,ch));
    end
end

% returns all ones, so is normally distributed

%% Calculate p value and correlation's z-stat against the null:
clear Pvalue zRho
for ss = 1:17
    for ch = 1:96
        clear NullDist Rho 
        NullDist = zvaluesScram.act_all_z(:,ss,ch);
        
        % All trials
        Rho = zvalues.allt_alltrials_z{ss,ch};
        
        % p value 
        Pvalue(ss,ch)=sum(abs(NullDist(:))>=abs(Rho))/numel(NullDist);% two-sided
        
        % z-statistic against the null distribution
        z = @(p,sign_z) sign_z .* norminv(1.0 - p./2.0 - eps, 0, 1);
        zRho(ss,ch) = z(Pvalue(ss,ch), sign(Rho-median(NullDist)));
        
        % Specified Trials
        clear Rho_spec
        Rho_spec = zvalues.allt_alltrials_z_spec{ss,ch};
        Pvalue_spec(ss,ch)=sum(abs(NullDist(:))>=abs(Rho_spec))/numel(NullDist);% two-sided
        zRho_spec(ss,ch) = z(Pvalue_spec(ss,ch), sign(Rho_spec-median(NullDist)));
        % Free Trials
        clear Rho_free
        Rho_free = zvalues.allt_alltrials_z_free{ss,ch};
        Pvalue_free(ss,ch)=sum(abs(NullDist(:))>=abs(Rho_free))/numel(NullDist);% two-sided
        zRho_free(ss,ch) = z(Pvalue_free(ss,ch), sign(Rho_free-median(NullDist)));
        
    end
end



%% Stats - NONPARAMETRIC
% define pttest
clear pttest
pttest.all_action_zero = zeros(5,ROInum);
%pttest.all_spec_free   = zeros(5,ROInum);

correction = 1;
for ch = 1:ROInum
    % ALL TIME
    
    % Action vs zero
    for ss = 1:length(sname); 
        action_allt(ss) = zRho(ss,ch);
    end
    [p,h,stats] = signrank(action_allt);
    pttest.all_action_zero(2,ch) = p;
    pttest.all_action_zero(1,ch) = h;
    pttest.all_action_zero(3,ch) = stats.zval;
    pttest.all_action_zero(4,ch) = mean(action_allt);
    pttest.all_action_zero(5,ch) = 0;
    
    
    % Specified vs Free
    for ss = 1:length(sname); spec_allt(ss) =  zRho_spec(ss,ch);
    free_allt(ss) =  zRho_free(ss,ch);end
    [p,h,stats] = signrank(spec_allt', free_allt', 'alpha', 0.05/correction, 'method', 'approximate');
    pttest.all_spec_free(2,ch) = p;
    pttest.all_spec_free(1,ch) = h;
        pttest.all_spec_free(3,ch) = stats.zval;
    pttest.all_spec_free(4,ch) = mean(spec_allt);
    pttest.all_spec_free(5,ch) = mean(free_allt);
    
    
    
end
 pval = pttest.all_action_zero(2,:)';
 pval_spec_free = pttest.all_spec_free(2,:)';

