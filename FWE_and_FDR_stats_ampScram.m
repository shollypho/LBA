
%% FDR
%% For 12-30Hz, Phase Scrambled

pttest_FDR = pttest;

pval_temp = ([0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;1;
1;0;0;0;0;0;1;1;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;1;
0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;1;0;1;0;0;1;1;1;1;0;0;0;0;0;0;0;0;
0;0;0;1;0;1;1;1;0;0;0;0;0;0;0;1]);

pttest_FDR.all_action_zero(1,:) = pval_temp;

% Spec vs free
pval_temp = zeros(96,1);
pttest_FDR.all_spec_free(1,:) = pval_temp;
pttest_FDR.all_spec_free_thresh = pttest.all_spec_free;

for ch = 1:96
pttest_FDR.all_spec_free_thresh(1,ch) = pttest_FDR.all_action_zero(1,ch)*pval_temp(ch,1);
end
pttest = pttest_FDR;

%% For 30-60Hz, AMP Scrambled

pttest_FDR = pttest;

pval_temp = ([0;1;1;1;0;1;1;1;1;0;0;1;1;1;0;1;
1;1;1;1;1;1;1;1;0;1;1;0;1;1;1;0;0;1;1;1;1;1;1;1;0;1;1;1;1;1;1;1;
0;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;0;0;1;1;1;1;
0;1;1;1;1;1;1;1;0;1;1;1;1;1;1;1]);

pttest_FDR.all_action_zero(1,:) = pval_temp;



% Spec vs free
pval_temp = zeros(96,1);
pttest_FDR.all_spec_free(1,:) = pval_temp;
pttest_FDR.all_spec_free_thresh = pttest.all_spec_free;

for ch = 1:96
pttest_FDR.all_spec_free_thresh(1,ch) = pttest_FDR.all_action_zero(1,ch)*pval_temp(ch,1);
end
pttest = pttest_FDR;

%% For 60-90Hz, Amp Scrambled

pttest_FDR = pttest;

pval_temp = ([1;1;1;1;1;1;1;1;0;0;0;0;1;0;0;0;
1;1;1;1;1;1;1;1;0;1;1;0;1;1;1;1;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;0;
0;1;1;1;1;1;1;1;1;1;1;1;1;0;1;0;1;1;1;1;1;1;0;1;0;1;1;1;1;1;1;1;
1;0;0;1;0;0;0;0;1;1;1;1;1;1;0;0]);

pttest_FDR.all_action_zero(1,:) = pval_temp;

% Spec vs free
pval_temp = zeros(96,1);
pttest_FDR.all_spec_free(1,:) = pval_temp;
pttest_FDR.all_spec_free_thresh = pttest.all_spec_free;

for ch = 1:96
pttest_FDR.all_spec_free_thresh(1,ch) = pttest_FDR.all_action_zero(1,ch)*pval_temp(ch,1);
end
pttest = pttest_FDR;
