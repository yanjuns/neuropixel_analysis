function [spatialcell_sec, spatialcell_spike, meanFRd, Tinfo] = spatial_info_shuffle(S, T, speed, post, posx, sp, trial, TrackEnd, trial_range, nboot, occThresh, p)
% function to shuffle spikes and calculate spatial information to identify
% spatial cells
% inputs
%     S, T
%     speed: animal's running speed (real coords), obtained from calcSpeed
%     post: time stamps
%     posx: position over time
%     sp: structure of data
%     trial
%     TrackEnd: track length
%     trial_range: specify baseline trials
%     p: params (time bin, etc)
% outputs:
%     spatialcell_sec: spatial cells that defined by bit/sec
%     spatialcell_spike: spatial cells that defined by bit/spike
%     Tinfo: a table contains above information
% Yanjun Sun, Stanford University, 2/8/2020
%% setup parameters
if ~exist('p', 'var') || isempty(p)
    SpeedCutoff = 2; %cm/s
    TrackStart = 0;
    SpatialBin = 5;
    TimeBin = median(diff(post));
else
    SpeedCutoff = p.SpeedCutoff;
    TrackStart = p.TrackStart;
    SpatialBin = p.SpatialBin;
    TimeBin = p.TimeBin;
end
if ~exist('occThresh', 'var') || isempty(occThresh)
    occThresh = 0.1;
end
if ~exist('nboot', 'var') || isempty(nboot)
    nboot = 1000;
end
if ~exist('trial_range', 'var') || isempty(trial_range)
    trial_range = [1:50];
end
if ~exist('TrackEnd', 'var') || isempty(TrackEnd);
    TrackEnd = 320; %John Wen's normal track
end
%% calculate normal spatial information
% trial_range = [1:100]; % trials to define, using which trials to define place cells
nCells = length(S);
Sd = {}; FRd = {};
for ii = 1:nCells
    Sd{ii} = S{ii}(trial_range,:);
end
Sd = cellfun(@sum, Sd, 'UniformOutput', false);
Td = sum(T(trial_range,:));
FRd = cellfun(@(x) x./Td, Sd, 'UniformOutput', false);
meanFRd = (cellfun(@(x) sum(sum(x)), Sd) / sum(sum(Td)))';
bitpsec = []; bitpspike = [];
for ii = 1:length(Sd);
    [bitpsec(ii,1), bitpspike(ii,1)] = spatial_info(FRd{ii}, meanFRd(ii), Td, occThresh);
end

%% prepare for filtering slow speed frames
%select trials
trialstart = find(trial == min(trial_range), 1);
trialend = find(trial == max(trial_range), 1, 'last');
post = post(trialstart:trialend);
posx = posx(trialstart:trialend);
speed = speed(trialstart:trialend);
%filter slow speed frames
% posxf = posx; posxf(speed <= SpeedCutoff) = NaN; 
slowf = find(speed <= SpeedCutoff);

%% to extract spiketime from all neurons
cells_to_plot = sp.cids(sp.cgs==2); 
spike_tall = {};
for jj = 1:nCells
    spike_t1 = sp.st(sp.clu==cells_to_plot(jj));
    spike_t1 = spike_t1(spike_t1 <= max(post) & spike_t1 >= min(post));
    [~,~,spidx] = histcounts(spike_t1,post);
    Lia = ~ismember(spidx, slowf);
    spike_tpass = spike_t1(Lia);
    spike_tall{jj} = spike_tpass;
end

%% shuffle
bitpsecboot = []; bitpspikeboot = [];
deltaT = randi([20,round(max(post)-min(post))-20],nboot,1);
binedges = TrackStart:SpatialBin:TrackEnd;
Tboot = repmat(Td,nCells,1);
tic
parfor iter = 1:nboot
    dT = deltaT(iter);
        Sboot = [];
        for k = 1:nCells
            spike_t = spike_tall{k};
            spike_tboot = spike_t + dT;
            idx = spike_tboot > max(post);
            spike_tboot(idx) = spike_tboot(idx) - max(post) + min(post);
            spike_t = sort(spike_tboot); % sort shuffled spikes
            [~,~,spike_idx] = histcounts(spike_t,post);
            % count spikes into bins
            spikeCount = histcounts(posx(spike_idx), binedges);
            Sboot(k,:) = spikeCount;
        end
        % cal shuffled spatial information
        FRboot = Sboot./Tboot;
        meanFRboot = sum(Sboot,2)./sum(Tboot,2);
        for jj = 1:nCells;
        [bitpsecboot(jj,iter), bitpspikeboot(jj,iter)] = spatial_info(FRboot(jj,:), meanFRboot(jj), Td, occThresh);
        end
end
toc
%% output final results as a table
bitpsecthresh = quantile(bitpsecboot, 0.99, 2);
bitpspikethresh = quantile(bitpspikeboot, 0.99, 2);
spatialcell_sec = find(bitpsec > bitpsecthresh);
spatialcell_spike = find(bitpspike > bitpspikethresh);
% output final results into a table
allneuron = [1:nCells]';
Lia1 = ismember(allneuron,spatialcell_sec);
Lia2 = ismember(allneuron,spatialcell_spike);
Tinfo = table(allneuron,cells_to_plot',bitpsec,Lia1,bitpspike,Lia2,'VariableNames',...
    {'neuron','cellID','bitpsec','spatialCell_sec','bitpspike','spatialCell_spike'});
% Tinfo = sortrows(Tinfo,{'bitpspike'},{'descend'});

end