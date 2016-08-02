function [X] = wmDrop_makeX1(stimLocs,stimSize,rfPtsX,rfPtsY,fwhm,res)
% TCS, 10/24/14 - implements build_basis_pts.m, make_stim_mask.m


% FOV will be defined by the basis set - for computing encoding model, at
% least
%
%
% note: does not!!! normalize X, need to do that wherever this is called
% (channelRespAmp)





% go from FWHM to rad
filt_size =  (1/rad2fwhm(1)) * fwhm;

% generate x, y grid for make_stim_mask and build_basis
% only need to go to the edge of the most distant stimulus x/y (as all
% other stim_mask values will be zero)
fov = 2*max(stimLocs(:)) + 2*stimSize;

% maybe add checks here to make sure spacing == fwhm

gridPts = linspace(-fov/2,fov/2,res);
[gridx, gridy] = meshgrid(gridPts,gridPts);

gridx = reshape(gridx,numel(gridx),1);
gridy = reshape(gridy,numel(gridy),1);
% build our basis set (nChannels x nPixels)
basis_set = build_basis_pts(rfPtsX,rfPtsY,filt_size,gridx,gridy);
% numel(gridx) x length(rfPtsX);



stim_mask = make_stim_mask(stimLocs(:,1),stimLocs(:,2),stimSize,gridx,gridy);

% now, we need to make a stimulus filter for each trial
% stim.xLocDeg, yLocDeg contains all center stim position in screen coords
% p.radDeg is size of stimulus

% stim_mask = nan(length(stim.xLocDeg),res(1)*res(2));
% 
% for ss = 1:length(stim.xLocDeg)
%     
%     thisr = sqrt((gridx-stim.xLocDeg(ss)).^2 + (gridy-stim.yLocDeg(ss)).^2);
%     thisstim = thisr <= p.radDeg;
%     
%     stim_mask(ss,:) = reshape(thisstim,1,numel(thisstim));
%     clear thisstim thisr;
% end

% stim_mask is nTrials x nPixels

% we want nTrials x nChannels

X = stim_mask' * basis_set;

% intentionally not normalizing X here, it should be normalized across all
% runs, and this will work on a single run...

return