function [Sindiv_block, Tindiv_block, FRindiv_block, FRSindiv_block] = calc_spatial_firingrate_block(Sindiv, Tindiv, trials_per_block, p)
% This function aims to calculate averaged firing rate in spatial bins 
% across multiple tirals. 
% Inputs are the ouputs from calc_spatial_firingrate:
% Sindiv: spike counts for each indivdual neuron
% Tindiv: occupancy time in each spatial in for each individual neuron
% Created Yanjun Sun, 8/8/2019

if ~exist('p', 'var') || isempty(p);
    SpatialBin = 5;
    SmoothSigmaFR = 15;
else
    SpatialBin = p.SpatialBin;
    SmoothSigmaFR = p.SmoothSigmaFR;
end

if ~exist('trials_per_block', 'var') || isempty(trials_per_block);
    trials_per_block = 10;
end
numTrial = size(Sindiv,1);
Sindiv_block = []; Tindiv_block = [];
for ii = 1:ceil(numTrial/trials_per_block);
    if ii < ceil(numTrial/trials_per_block);
        Sindiv_block(ii,:) = nansum(Sindiv(trials_per_block*ii - (trials_per_block-1):trials_per_block*ii, :));
        Tindiv_block(ii,:) = nansum(Tindiv(trials_per_block*ii - (trials_per_block-1):trials_per_block*ii, :));
    else
        Sindiv_block(ii,:) = nansum(Sindiv(trials_per_block*ii - (trials_per_block-1):end, :));
        Tindiv_block(ii,:) = nansum(Tindiv(trials_per_block*ii - (trials_per_block-1):end, :));
    end
end
FRindiv_block = Sindiv_block./ Tindiv_block;

% gaussian filter for smoothing
smoothSigma = SmoothSigmaFR/SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

% smooth firing rate
FRSindiv_block = [];
for jj = 1:size(FRindiv_block,1);
    firing_rate_smooth = conv(repmat(FRindiv_block(jj,:),1,3),gauss_filter,'same');
    FRSindiv_block(jj,:) = firing_rate_smooth(numel(FRindiv_block(jj,:))+1:numel(FRindiv_block(jj,:))*2);
end
end
