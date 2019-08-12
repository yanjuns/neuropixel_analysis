function [spikeCount, occTime, firing_rate, firing_rate_sth] = calc_spatial_firingrate(idx, posx, post, trial, ntrial, p, TrackEnd)
% calculates smoothed firing rate on linear track
% Malcolm Campbell 5/21/15
% modified 6/6/18 MGC
% all time bins are now equal length
% edited kei masuda 7/1/19
% edited Yanjun Sun 8/7/19 with multiple outputs
% inputs:
%     idx: path to neuropixel .mat file
%     p: params (spatial bin size, etc)
%     TrackEnd: typically, the length of the track 
% outputs:
%     spikeCount: number of spikes in each spatial bin
%     occTime: occupancy time for each spatial bin
%     firing_rate: firing rate over position
%     firing_rate_sth: smoothed firing rate
%     divide spike counts by occupancy

if ~exist('TrackEnd', 'var') || isempty(TrackEnd);
    TrackEnd = 640; %John Wen's extened track
end

if ~exist('p', 'var') || isempty(p);
    TrackStart = 0;
    SpatialBin = 5;
    TimeBin = median(diff(post));
    SmoothSigmaFR = 15;
else
    TrackStart = p.TrackStart;
    SpatialBin = p.SpatialBin;
    TimeBin = p.TimeBin;
    SmoothSigmaFR = p.SmoothSigmaFR;
end

binedges = TrackStart:SpatialBin:TrackEnd;
time_per_bin = histcounts(posx(trial == ntrial), binedges);
occTime = time_per_bin * TimeBin;
spikeCount = histcounts(posx(idx), binedges);
firing_rate = spikeCount./occTime;

% % interpolate missing values
% if sum(isnan(firing_rate))>0
%     firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
% end

% gaussian filter for smoothing
smoothSigma = SmoothSigmaFR/SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

% smooth firing rate
firing_rate_smooth = conv(repmat(firing_rate,1,3),gauss_filter,'same');
firing_rate_sth = firing_rate_smooth(numel(firing_rate)+1:numel(firing_rate)*2);

end
