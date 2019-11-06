function [fChan,fChanRS] = envelope_ROI(data, fHpi, fLpi)

%combine planar gradiometers
%D.chantype gives the type of sensors: they are always distributed in
%triplets of the type {Megplanar(X-coord) Megplanar(Y-coord) MEGMag}
%basically the  Megplanar x/y and megmag sensors belonging to a triplet lie on the
%same point in the scalp.
%The traditional way to combine x/y Megplanar into a scalar is by taking
%their RMS. However it poses some problems when conditions differ in number of
%trials. An alternative is to average them.

 
 %NOTE:Data are stored in a 3-D matrix: channels x samples x trials;
 % average pairs of Megplanar channels
 
%
%% Filter data in the desired range and calculate Hilbert envelope
clear fChan fChanRS

D.fsample = 250;

for ii = 1:size(data,3)
    
    fChan(:,:,ii) = ft_preproc_bandpassfilter(data(:,:,ii), D.fsample, [fHpi fLpi], 4, 'but', 'twopass');
    fChan(:,:,ii) = abs(hilbert(fChan(:,:,ii)'))';
end

% Combine gradiometers by using the mean: note this makes only sense for
% the power which is always positive.
% if strcmp(ChanType,'MEGPLANAR');
%     fChanRS = squeeze(mean(reshape(fChan,2,102,size(fChan,2),size(fChan,3))));
%     fChanRS_tr = fChanRS;
%     fChanRS = mean(fChanRS, 3);
%     fChanmean = mean(fChan,3);
% else
    fChanRS = mean(fChan,3);
% end


