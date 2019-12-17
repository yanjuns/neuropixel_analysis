function [S, T, FR, FRS, corrMatrix, S_block, T_block, FR_block, FRS_block, corrBlock, cells_to_plot, speed]...
    = get_ratemap_n_corr(matPath, trackLength, trials_per_block, filterspeed, speedthresh, paramsPath)

% John Wen 7/1/19
% Kei Masuda 7/3/19
% Yanjun Sun 8/8/2019
% Runs Malcolm's firing rate code on all good cells within a .mat file. 
% Calculate firing rate by dividing spike counts by occupancy
% inputs:
%     matPath: path to .mat file after running sync_vr_to_np.m. Specify as
%              string.
%              e.g. '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190620_keicontrasttrack_ketamine1_g0/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine'
%     trackLength: specify the length of the track
%     trials_per_block: specify average from how many trials
%     paramsPath: Optional argument. Path to parameters file, which includes 
%                 spatial bin size and smoothing kernel. 
%                 Default = '/Volumes/groups/giocomo/export/data/Projects/ ...
%                 JohnKei_NPH3/UniversalParams'
% outputs:
%     S: spike counts for each cell
%     T: occupancy time for the mouse
%     FR: firing rate for each cell
%     FRS: smoothed firing rate for each cell
%     corrMatrix: all the correlation matrix for each individual trial each
%     neuron

if ~exist('speedthresh', 'var') || isempty(speedthresh)
    speedthresh = 2; %cm/s
end
if ~exist('filterspeed', 'var') || isempty(filterspeed)
    filterspeed = true;
end

%% Load .mat and params files
% load specific components from .mat file
load(fullfile(matPath), 'post','posx','sp','trial'); 

% load params file
if ~exist('paramsPath','var')
%     paramsPath = '\\oak-smb-giocomo.stanford.edu\groups\giocomo\ysun\Neuropixels\UniversalParams.xlsx';
    params = readtable('UniversalParams.xlsx'); %make sure the file is on the path
else
    params = readtable(paramsPath);
end

%% calculate firing rates for only good cells
cells_to_plot = sp.cids(sp.cgs==2); 
nCells = size(cells_to_plot, 2);

%% Calculate firing rates
% specify inputs into the firing rates calculation
trackEnd = trackLength;
p = params;

% calculate the firing rate for a single cell across all trials
S = cell(1,nCells);
T = [];
FR = cell(1,nCells);
FRS = cell(1,nCells);
corrMatrix = [];
S_block = cell(1,nCells);
T_block = [];
FR_block = cell(1,nCells);
FRS_block = cell(1,nCells);
corrBlock = [];
numTrial = max(trial);
speed = calcSpeed(posx,post);
%prepare for filtering slow speed frames
posxf = posx; posxf(speed <= speedthresh) = NaN; 
slowf = find(speed <= speedthresh);

for k = 1:nCells
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
    % get spike times and index into post for cell k 
    spike_t = sp.st(sp.clu==cells_to_plot(k)); 
    [~,~,spike_idx] = histcounts(spike_t,post);
%     [~,~,spikeTrial_idx] = histcounts(spike_t,trial); % NOT SURE IF WE STILL NEED THIS
    % filter out slow speed frame if needed
    if filterspeed
        Lia = ismember(spike_idx, slowf);
        spike_idx(Lia) = [];
        % for cell k, iteratively calculate the firing rate for each trial
        Sindiv = []; Tindiv = []; FRindiv = []; FRSindiv = [];
        for ii = 1:numTrial
            idx = spike_idx(trial(spike_idx)==ii);
            [spikeCount, occTime, firing_rate, firing_rate_sth]...
                = calc_spatial_firingrate(idx, posxf, post, trial, ii, [], trackEnd);
            Sindiv(ii,:) = spikeCount;
            Tindiv(ii,:) = occTime;
            FRindiv(ii,:) = firing_rate;
            FRSindiv(ii,:) = firing_rate_sth;
        end
        S{k} = Sindiv;
        T = Tindiv;
        FR{k} = FRindiv;
        FRS{k} = FRSindiv;
        % calculate trial by trial correlations for one cell 
        corrMatrix(:,:,k) = corr(FRSindiv');
        % calculate correlations across blocks of trials (every 10 trials)
        [Sindiv_block, Tindiv_block, FRindiv_block, FRSindiv_block]...
            = calc_spatial_firingrate_block(Sindiv, Tindiv, trials_per_block);
        S_block{k} = Sindiv_block;
        T_block = Tindiv_block;
        FR_block{k} = FRindiv_block;
        FRS_block{k} = FRSindiv_block;
        % calculate trial block by trial block correlations for one cell 
        corrBlock(:,:,k) = corr(FRSindiv_block');
    else
        Sindiv = []; Tindiv = []; FRindiv = []; FRSindiv = [];
        for ii = 1:numTrial
            idx = spike_idx(trial(spike_idx)==ii);
            [spikeCount, occTime, firing_rate, firing_rate_sth]...
                = calc_spatial_firingrate(idx, posx, post, trial, ii, [], trackEnd);
            Sindiv(ii,:) = spikeCount;
            Tindiv(ii,:) = occTime;
            FRindiv(ii,:) = firing_rate;
            FRSindiv(ii,:) = firing_rate_sth;
        end
        S{k} = Sindiv;
        T = Tindiv;
        FR{k} = FRindiv;
        FRS{k} = FRSindiv;
        % calculate trial by trial correlations for one cell 
        corrMatrix(:,:,k) = corr(FRSindiv');
        % calculate correlations across blocks of trials (every 10 trials)
        [Sindiv_block, Tindiv_block, FRindiv_block, FRSindiv_block]...
            = calc_spatial_firingrate_block(Sindiv, Tindiv, trials_per_block);
        S_block{k} = Sindiv_block;
        T_block = Tindiv_block;
        FR_block{k} = FRindiv_block;
        FRS_block{k} = FRSindiv_block;
        % calculate trial block by trial block correlations for one cell 
        corrBlock(:,:,k) = corr(FRSindiv_block');
    end
end

%plot for averaged correlation
avg_corrMatrix = nanmean(corrMatrix,3);
avg_corrblock = nanmean(corrBlock,3);
figure;
imagesc(avg_corrMatrix);
axis square;
% colormap pink;
xlabel('Trials');
saveas(gcf,'singleTrialCorr.fig');

figure;
imagesc(avg_corrblock);
axis square;
% colormap pink;
xlabel('Trials chunks');
saveas(gcf,'trialChunkCorr.fig');

end