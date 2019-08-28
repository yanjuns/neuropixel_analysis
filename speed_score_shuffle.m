function speedthresh = speed_score_shuffle(speed, post, sp, trial, trial_range, p)
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
nboot = 100;
deltaT = randi([30,round(max(post))-30],nboot,1);
cells_to_plot = sp.cids(sp.cgs==2); 
nCells = size(cells_to_plot, 2);
speedshuff = [];
parfor iter = 1:nboot
    speedboot = [];
    dT = deltaT(iter);
    for k = 1:nCells
        spike_t = sp.st(sp.clu==cells_to_plot(k));
        spike_t = spike_t(spike_t <= max(post));
        spike_tboot = spike_t + dT;
        idx = spike_tboot > max(post);
        spike_tboot(idx) = spike_tboot(idx) - max(post);
        spike_t = sort(spike_tboot); % sort shuffled spikes
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
        speedboot = [speedboot; speedScore];
    end
    speedshuff = [speedshuff;speedboot];
end
%% output final results as a table
speedthresh = [];
speedthresh(1,1) = quantile(speedshuff, 0.01);
speedthresh(1,2) = quantile(speedshuff, 0.99);
figure;
histogram(speedshuff);
hold on;
line([speedthresh(1,2), speedthresh(1,2)], get(gca,'ylim'), 'LineStyle','--','Color','r','LineWidth',1);
text(speedthresh(1,2)+0.005,range(ylim)/1.5,['Thresh = ',num2str(speedthresh(1,2),'%.3f')],'Color','red','FontSize',8)
line([speedthresh(1,1), speedthresh(1,1)], get(gca,'ylim'), 'LineStyle','--','Color','g','LineWidth',1);
text(speedthresh(1,1)+0.005,range(ylim)/1.5,['Thresh = ',num2str(speedthresh(1,1),'%.3f')],'Color','green','FontSize',8)

end