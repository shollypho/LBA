% Scramble phase for each trial and record r values
% Permute 1000 times
% clear all % be careful with clear all b/c trialdata takes a long time to
% load
%% Set up
%addpath(genpath('/imaging/hp02/spm12b'));
addpath(genpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/ERPs/average_variable_ndt'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt');
addpath('/imaging/hp02/software_n_scripts');

sname = [23 24 25 26 27 28 29 30 31 32 33 527 528 529 530 533 534];

ROInum=96;

% load('trialdata.mat'); % commented out for now b/c it takes a long time to load
% load('lba_stats.mat');
%load('ndt_split.mat');
%load('rvalues');
disp('check1')
%% open matlabpool if required
% %     ParType = 0;  % Fun on Login machines (not generally advised!)
% %     ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
ParType = 2;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)


if ParType
    if matlabpool('size')==0;
        %MaxNsubs = 1;
        %if ParType == 2
        %    for g=1:length(cbu_codes)
        %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
        %    end
        %end
        P = cbupool(96);
        P.ResourceTemplate='-l nodes=^N^,mem=16GB,walltime=72:00:00';
        matlabpool(P);
    end
end

%% Extract timeseries data from each trial/channel/subject, scramble and correlate with LBA
% Need to do for both all time and acc
disp('check2')
%set up rvalueScram cell
% with LBA
full_r_ch = zeros(10000,length(sname), ROInum); 
full_p_ch = zeros(10000,length(sname), ROInum);
             
parfor ch = 1:ROInum
    ch
    all_r_mean = zeros(length(sname),1);
    all_p_mean = zeros(length(sname),1);
    full_r = zeros(10000,length(sname)); full_p = zeros(10000,length(sname));
 
    for ss = 1:length(sname)
        ss
        % Grab ndt split for this subject and channel
        ndt = floor(trialdata{ss,ch}.ndt_samp);
        t1 = round(all_ndt_split(ss, ch));
        if t1 ==0; t1 = 1;end
        if t1>ndt; t1 = ndt; end
        
        tr_num = length(trialdata{ss,1}.trial);
        
        cnt = 0;
        
        % concatinate all trials
        
        
        fullmodel = []; full_data = [];
        for tr = 1:length(trialdata{ss,1}.trial)
            
            
            % LBA timeseries (keep in samples for ease)
            if ~isempty(trialdata{ss,ch}.trial{1,tr})
                RT = trialdata{ss,ch}.trial{1,tr}.RT_samp;
                b1 = trialdata{ss,ch}.trial{1,tr}.LBA_grad;
                t2 = RT - (ndt-t1);
                acc_time = RT-ndt;
                thres = b1*acc_time;
                % original timeseries which has been correlated with LBA,
    
                fullmodel = [fullmodel zeros(1,t1-1) 0:b1:thres ones(1,RT-t2)*(thres)]; % model from t=0 to t = RT
                full_data = [full_data trialdata{ss,ch}.trial{1,tr}.trialdata];
            end
        end
      
    
        
        % Fast fourier transform to get phase
        fft_full = fft(full_data);
        f1 = abs(fft_full);
        f2 = angle(fft_full); % or imag?? PHASE
        
        
        % permute through scramble of phase, inverse FFT and correlate
        
        for perm_phase = 1:10000
            % scramble the phase
            
            f2r = randswap(f2);
            fft_full_scram = f1.*exp(1i*f2r);

            % Invert back to time domain
            ifft_full_scram = abs(real(ifft(fft_full_scram)));
            
            % correlate with LBA model
            
            [rho,pval] = corr(fullmodel', ifft_full_scram', 'Type', 'Spearman');

            full_r(perm_phase,ss) = rho; full_p(perm_phase,ss) = pval;
            
        end
        
    end
    full_r_ch(:,:,ch) = full_r;
    full_p_ch(:,:,ch) = full_p;
    
end
% Put into one struct.... parfor wouldn't let me do this with the parfor

mean_full_r =squeeze(mean(full_r_ch,1));
mean_full_p =squeeze(mean(full_p_ch,1)); 

rvaluesScram.all_r_mean = mean_full_r;
rvaluesScram.all_p_mean = mean_full_p;
rvaluesScram.all_r = full_r_ch;
rvaluesScram.all_p = full_p_ch;


save('60-90Hz/rvaluesScram_phase_spearman', 'rvaluesScram');

%close cluster
% if ParType
%     matlabpool close
% end