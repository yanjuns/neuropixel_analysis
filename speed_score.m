function speedinfo = speed_score(speed, post, sp, trial, trial_range, p)
% function to calculate speed score in VR using linear fit
% inputs
%     speed: animal's running speed (real coords), obtained from calcSpeed
%     post: time stamps
%     sp: structure of data
%     trial
%     trial_range: specify baseline trials
%     p: params (time bin, etc)
% outputs:
%     speedScore: speed score
%     params: [intercept, slope] of linear fit
%     speedinfo: a table contains above information
%
% based on Kropff et al 2015
% Malcolm Campbell, 3/30/16, modified 6/6/18
% Yanjun Sun, 8/26/2019
%% setup parameters
if ~exist('p', 'var') || isempty(p)
    SpeedCutoff = 2; %cm/s
    TimeBin = median(diff(post));
else
    SpeedCutoff = p.SpeedCutoff;
    TimeBin = p.TimeBin;
end
smoothSigma = 20;

%% calculate instantanenous speed in specified trials
trialstart = find(trial == min(trial_range), 1);
trialend = find(trial == max(trial_range), 1, 'last');
post = post(trialstart:trialend);
speed = speed(trialstart:trialend);
select = speed > SpeedCutoff; 
% speed = reshape(speed,numel(post),1);
speedFilt = speed(select);

%% calculate instantaneous fr in specified trials for each neuron
cells_to_plot = sp.cids(sp.cgs==2); 
nCells = size(cells_to_plot, 2);
speedinfo = [];
for k = 1:nCells
    spike_t = sp.st(sp.clu==cells_to_plot(k)); 
    % estimate instantaneous firing rate by smoothing spike count histogram
    h = histcounts(spike_t,post);
    h = [0,h];
    fr = reshape(h,numel(post),1);
    fr = fr./TimeBin;
    fr = gauss_smoothing(fr,smoothSigma);
    % throw out bins with speeds below threshold
    frFilt = fr(select);
    % compute speed score by correlating instantaneous firing rate and speed
    speedScore = corr(speedFilt,frFilt);
    params = [ones(sum(select),1) speedFilt]\frFilt;
    speedinfoindiv = [cells_to_plot(1,k), speedScore, params'];
    speedinfo = [speedinfo; speedinfoindiv];
end

%% output final results as a table
speedinfo = array2table(speedinfo,'VariableNames',{'CellID','SpeedScore','Intercept','Slop'});

end