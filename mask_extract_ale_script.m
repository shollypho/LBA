%% Extracts timecourses from vertices within an ROI mask image
% Alex C (10/2012)
% Modificat,ions Ece K (10/2015)
% Modifications Holly P (02/2018)
% Change to include Alessandro Tomassini's method

% General variable setup
%addpath('/work/imaging8/MEG/bananas_2015_EK/scripts/NIfTI_20140122/');
addpath /imaging/local/meg_misc/
addpath /hpc-software/matlab/cbu/;
addpath /imaging/hp02/software_n_scripts/NIfTI_20140122
addpath /imaging/hp02/ece_scripts

%clear all; clc;
datadir = '/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/parcellations_AleT_100Hzds_m500_1500ms/';
maskdir='/imaging/hp02/finger_tapping08/analysis_spm/LBA_modelling/Time-Frequency/average_variable_ndt/mat_files/source_reconstruction/Harvard_Oxford/';

subs = {'s23','s24','s25','s26','s27','s28','s29','s30','s31','s32','s33',...
    's527','s528','s529','s530','s533','s534'};

HO_list; % list of Harvard Oxford parcellations

% Filename of ROI masks and nicname used in naming new files
maskname = HO_list_long;
masknic = HOList_nickname; %

homo = 0; % do RH homologues of LH ROIs?
reslice=1;
newvoxel=[1 1 1];

% Parts of filename with source localised info in
front = 'devracbefffM';
lastbit = '_FT_blk1.mat';

% Options
wholebrain=0; % normal mesh 8196 vertices. if you want more you need to increase your mesh size
usemask =1;
virtualelectrode=0;% use ROI or single coord
val = 1;  % inversion number
virtual_el_coor=[];

%% Get coordinates of vertices in ROIs

%% open matlabpool if required
%ParType = 0;  % Fun on Login machines (not generally advised!)
%ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
ParType =0;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)


if ParType
    if matlabpool('size')==0;
        %MaxNsubs = 1;
        %if ParType == 2
        %    for g=1:length(cbu_codes)
        %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
        %    end
        %end
        P = cbupool(17);
        matlabpool(P);
    end
end
%%

newsample = 100;%Hz

for ss = 1:length(subs)
    %% ROIs identification
    
    % load source reconstructed data
    D = spm_eeg_load(sprintf('%s/%s/%s%s%s', datadir, subs{ss}, front, subs{ss}, lastbit));
    
    NII = load_nii('/imaging/at07/templates/HarvardOxford-cort-maxprob-thr25-2mm.nii');
    
    T = [NII.hdr.hist.srow_x; NII.hdr.hist.srow_y; NII.hdr.hist.srow_z; [0 0 0 1]];
    
    % convert MNI to Voxels
    mni2vox = @(xyz) [xyz 1] * inv(T)' + 1;
    
    idmin = [];
    
    for iv = 1:size(D.inv{invnum}.mesh.tess_mni.vert,1)
        
        vox = round(mni2vox(D.inv{invnum}.mesh.tess_mni.vert(iv,:)));
        idmin(iv,1) = NII.img(vox(1),vox(2),vox(3));%#ok
        %idmin(iv,2) = (vox(1)<45);%#ok 0Left 1Right :split ROI in two lateral halves to create binarty ROIs
        idmin(iv,2) = D.inv{invnum}.mesh.tess_mni.vert(iv,1)>0;%#ok0Left 1Right
    end
    
    %compute inversion Weights
    ind = {};
    for m=1%:3
        ind{m} = 1:size(D.inv{invnum}.inverse.U{m},1);%#ok
        if m>1
            ind{m} = ind{m}+ind{m-1}(end);%#ok
        end
    end
    
    %W = {};
    W=[];
    
    for m=1%:3
        W{m} = D.inv{invnum}.inverse.M(:,ind{m}) * D.inv{invnum}.inverse.U{m};%#ok
        %    I = [I D.inv{1}.inverse.M(:,ind{m}) * D.inv{1}.inverse.U{m}];
    end
    
    
    W = [W{1}];% W{2} W{3}];%weight matrix
    
    % select only the good channels
    indch = [D.indchantype('MEGPLANAR')];% D.indchantype('MEG') D.indchantype('EEG')];
    
    indch = setdiff(indch,D.badchannels);
    
    data = D(indch,:,:);
    
    %project sensors in source space
    %for i = 1:size(W,1)
    
    %w = W(i,:);
    newd = [];
    for itrial = 1:size(data,3)
        newd(:,:,itrial) = W*data(:,:,itrial);%#ok [1 x Nchans] * [Nchans x time]
    end
    
    %end
    
    %% Extract ROIs PCA
    %The prob thresh Harvard-Oxford atlas has indices from 0 (no match) to 48
    origs = size(newd);
    
    ROIlabs = unique(idmin(:,1));ROIlabs = ROIlabs(ROIlabs>0);
    newdata = []; varexp = [];
    for hem = 1:2
        
        for ROI = 1:numel(ROIlabs)
            ss
            ROI
            
            idROI = find(idmin(:,1)==ROIlabs(ROI) & idmin(:,2)==(hem-1));
            [~,sc,~,~,vxp] = pca(reshape(newd(idROI,1:201,1:length(D.condlist)-1),numel(idROI),[],1)','NumComponents',1);
            sc = sc';
            varexp(ROIlabs(ROI)) = vxp(1);%#ok
            newdata(ROIlabs(ROI),:,:,hem) = reshape(sc,1,201, length(D.condlist)-1);%origs(2),origs(3));%#ok
            
        end
    end
    
    newd_cat = [];
    newd_cat = cat(1,newdata(:,:,:,1),newdata(:,:,:,2));

    mask_1stEigen(ss).data = newd_cat;
    
    mask_1stEigen(ss).dInfo.time    = D.time;
    %dInfo.behav   = D.behav;
    mask_1stEigen(ss).dInfo.fsample = D.fsample;
    mask_1stEigen(ss).dInfo.condlist= D.condlist;
    mask_1stEigen(ss).dInfo.conds   = D.conditions;
    mask_1stEigen(ss).dInfo.ROIs    = idmin;
    mask_1stEigen(ss).dInfo.PCAvExp = varexp;
    mask_1stEigen(ss).dInfo.NII     = NII.fileprefix;
    mask_1stEigen(ss).dInfo.invN    = invnum;
    
    
end

%save('mask_1stEigen.mat', 'mask_1stEigen');

% for mask = 1:10%length(masknic)
%
%     mask
%     disp(['Working on mask: ' masknic{mask}]);
%     if usemask
%
%         if reslice==0
%             P = ([maskdir  maskname{mask}]);
%             [M XYZ] = spm_read_vols(spm_vol(P));
%             XYZ_mask = XYZ(:,find(M))'; clear P M XYZ V
%         elseif reslice==1
%             nii_fname=[maskdir  maskname{mask}];
%             reslice_nii(nii_fname,[maskdir 'r' num2str(newvoxel(1)) 'mm_' maskname{mask}], newvoxel);
%             [M XYZ] = spm_read_vols(spm_vol([maskdir 'r' num2str(newvoxel(1)) 'mm_' maskname{mask}]));
%             XYZ_mask = XYZ(:,find(M))'; clear P M XYZ V
%         end
%
%         if homo
%             x = XYZ_mask(:,1).*-1;
%             XYZ_mask(:,1) = x; clear x
%         end
%
%         load('/imaging/hp02/ece_scripts/MEG_source_vertices_MNI.mat');
%         vertices = int16(vertices);
%
%         % Get mask coords for vertices only
%         member = ismember(vertices,XYZ_mask, 'rows');
%         XYZ = double(vertices(member,:)); clear member vertices XYZ_mask
%         disp(['Number of vertices: ' num2str(size(XYZ,1))]);
%
%     elseif virtualelectrode
%         XYZ = virtual_el_coor;
%     elseif wholebrain
%         load('/imaging/hp02/ece_scripts/MEG_source_vertices_MNI.mat');
%         vertices = double(vertices);
%         XYZ=vertices;
%     end
%
%
%     for ss = 1%:length(subs)
%
%         disp(['Extracting data from subject ' subs{ss} ' ...']);
%         theData = [datadir '/' subs{ss} '/' front subs{ss} lastbit];
%
%         D = spm_eeg_load(theData);
%
%         D.val = val;
%         D.inv{val}.source.fname = [ datadir '/' subs{ss} '/tempROI' '_' front subs{ss} lastbit];
%         D.inv{val}.source.type  = 'evoked'; %'trials' % was originally trials, but crashed
%         D.inv{val}.source.rad = 0;
%         D.inv{val}.source.XYZ = XYZ;
%
% %         numClusters = size(XYZ,1);
% %         for r = 1:numClusters
% %             D.inv{val}.source.label{r} = sprintf('ROI %d %d %d chan%d',XYZ(r,1), XYZ(r,2), XYZ(r,3), r);  % Insert your names if you want
% %         end
%
%         [Ds, D] = spm_eeg_inv_extract(D);
%
%         % extract first eigenvariate
%         for tr = 1:size(Ds,3)
%
%             y=0; u =0; s =0; v =0; Dy =0;
%
%             y=squeeze(Ds(:,:,tr)); % one trial at a time b/c svd won't work with 3 or more dimensions
%             %y=reshape(y,[size(Ds,1),tlength*nepoch]);
%
%             [u, s, v] = svd(y*y'); % gets all the eiganvariates
%             u       = u(:,1); % take the first eiganvariate, explains most of the variance
%             Dy = y'*u;
%
%             mask_1stEigen(mask,ss).trials(tr,:)=Dy';
%
%         end
%
%     end
% end
%
%
% %save('mask_1stEigen.mat', 'mask_1stEigen')
%
