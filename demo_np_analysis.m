%% Neuropixel data analysis streamline
%% reformat data if necessary, OPTIONAL
matPath2reform = 'D:\Giocomo_Neuropixels_analysis\E2_190619_MEC\E2_190619_johnrepeatingtrack_train100+meth21.mat'
load(matPath2reform);
savepath = 'D:\Giocomo_Neuropixels_analysis\E2_190619_MEC_reformat'
[posx,post,trial,sp] = reformat_rawdata (savepath, posx, post, trial, sp, [11,60,101,120]);

%% define trials for each animal
bslTrials = [1:100]; %specify for each particular mouse 
methTrials = [101:200]; %specify for each particular mouse 

%% Get firing rate map and rate map correlation for each neuron
matPath = 'D:\Giocomo_Neuropixels_analysis\E1_190619_MEC_reformat\E1_190619_MEC_reformat.mat '
load(matPath);
trackLength = 320 %John Wen's extended track
trials_per_block = 10;
[S, T, FR, FRS, corrMatrix, S_block, T_block, FR_block, FRS_block, corrBlock, cellID, speed]...
    = get_ratemap_n_corr(matPath, trackLength, trials_per_block, true);
% get mean firing rate for each cell and each trial
meanFR = NaN(size(S,2), size(T,1));
for ii = 1:size(S,2);
    for jj = 1:size(T,1);
        meanFR(ii,jj) = sum(S{ii}(jj,:))/sum(T(jj,:));
    end
end
save('ratemap_n_corr.mat', 'S', 'T', 'FR', 'FRS', 'corrMatrix',...
    'S_block', 'T_block', 'FR_block', 'FRS_block', 'corrBlock', 'cellID', 'speed', 'meanFR', '-v7.3');
%% plot ratemap figures
plot_ratemap(FRS, cellID)
%print figure into eps if needed
figname = 'Cell1to30_ratemap_auto';
fig = openfig(figname);
print (fig, '-painters', '-depsc', [figname, '.eps']);
%plot a specific one if needed
figure;
cellnum = find(cellID == 171);
pcolor(FRS{cellnum});
shading flat
colormap pink
figure;
imagesc(corrBlock(:,:,cellnum));
axis image
%% plot mean firing rate and running speed
% calculate speed by trials
speedbytrial = NaN(1, max(trial));
speed(speed < 2) = NaN;
for ii = 1:max(trial)
    speedbytrial(:,ii) = nanmean(speed(trial == ii));
end
save('ratemap_n_corr.mat', 'speedbytrial', '-append');
% plot
figure;
plot(mean(meanFR));
xlabel('Trials');
ylabel('Avg Firing Rate of all neurons (Hz)');
hold on;
yyaxis right;
plot(speedbytrial);
ylabel('Speed (cm/s)');
saveas(gcf, 'meanFRnSpeed.fig');

% test if mean firing rate significantly changed
bls = mean(meanFR(:, bslTrials),2);
meth = mean(meanFR(:, methTrials),2);
% bls = mean(meanFR(nonspeedcell, bslTrials),2);
% meth = mean(meanFR(nonspeedcell, methTrials),2);
[p,h] = signrank(bls, meth)
% test speed difference
[h,p] = ttest2(speedbytrial(bslTrials), speedbytrial(methTrials))

%% look at the correlation between firing rate and speed
% % overall correlation
% corr(speedbytrial(:,bslTrials)', mean(meanFR(:,bslTrials))')
% corr(speedbytrial(:,methTrials)', mean(meanFR(:,methTrials))')
% % bin by bin correlation
% FR_speed_corr = [];
% binsize = 10; k = 0;
% for ii = 1:binsize:max(trial)
%     k = k+1;
%     if ii+binsize-1 < max(trial)
%         FR_speed_corr(k,:) = corr(speedbytrial(:,ii:ii+9)', mean(meanFR(:,ii:ii+9))');
%     else
%         FR_speed_corr(k,:) = corr(speedbytrial(:,ii:end)', mean(meanFR(:,ii:end))');
%     end
% end
% [h,p] = ttest2(FR_speed_corr(1:ceil(max(bslTrials)/10),:), FR_speed_corr(11:ceil(max(methTrials)/10),:))

%% analyze lick behavior
load(matPath);
[lick_idx, lickMat, lickAccByTrial] = plot_lick_behav(lickt, lickx, post, posx, trial);

%% identify speed cells
% calculate speedscore, intercept, and slop of speed related coding
speedinfo = speed_score(speed, post, sp, trial, bslTrials);
% shuffling spikes to get a threshold for determine speed cells
tic;
speedthresh = speed_score_shuffle(speed, post, sp, trial, bslTrials);
toc;
% identify speed cells
speedidx = (speedinfo.SpeedScore > speedthresh(1,2));
negspeedidx = (speedinfo.SpeedScore < speedthresh(1,1));
speedinfo.speedidx = speedidx;
speedinfo.negspeedidx = negspeedidx;
nonspeedcell = ~speedidx & ~negspeedidx;
save('speedcell.mat','speedinfo','speedidx','negspeedidx', 'nonspeedcell');

% calculate speedscore, intercept, and slop of speed related coding for
% meth trials
tic
speedinfo_meth = speed_score(speed, post, sp, trial, methTrials);
% shuffling spikes to get a threshold for determine speed cells
tic;
speedthresh_meth = speed_score_shuffle(speed, post, sp, trial, methTrials);
toc;
% identify speed cells
speedidx = (speedinfo_meth.SpeedScore > speedthresh_meth(1,2));
negspeedidx = (speedinfo_meth.SpeedScore < speedthresh_meth(1,1));
speedinfo_meth.speedidx = speedidx;
speedinfo_meth.negspeedidx = negspeedidx;
% nonspeedcell = ~speedidx & ~negspeedidx;
save('speedcell.mat','speedinfo_meth', 'speedthresh_meth', '-append');
toc
%% find speed and non-speed cells and plot meanFR
%plot meanFR vs speed
figure;
hold on;
subplot(2,2,1);
plot(mean(meanFR));
xlabel('Trials');
ylabel('Avg Firing Rate of all neurons (Hz)');
yyaxis right;
plot(speedbytrial);
ylabel('Speed (cm/s)');
title('All Neurons');
subplot(2,2,2);
plot(mean(meanFR(nonspeedcell,:)));
xlabel('Trials');
ylabel('Avg Firing Rate of nonspeed cells (Hz)');
yyaxis right;
plot(speedbytrial);
ylabel('Speed (cm/s)');
title('Non-speed Cells');
subplot(2,2,3);
plot(mean(meanFR(speedidx,:)));
xlabel('Trials');
ylabel('Avg Firing Rate of speed cells (Hz)');
yyaxis right;
plot(speedbytrial);
ylabel('Speed (cm/s)');
title('Speed Cells');
subplot(2,2,4);
plot(mean(meanFR(negspeedidx,:)));
xlabel('Trials');
ylabel('Avg Firing Rate of neg speed cells (Hz)');
yyaxis right;
plot(speedbytrial);
ylabel('Speed (cm/s)');
title('Negative Speed Cells');
saveas(gcf, 'meanFRnSpeedmod.fig');

% calculate averaged mean firing rate of selected baseline trials and meth trials
% load('ratemap_n_corr.mat')
% load('speedcell.mat')
% bslT = 41:90; % selected number of sub trials if needed
% methT = 101:120;
% % for all the neurons
% bslFR = mean(meanFR(:, bslT),2);
% methFR = mean(meanFR(:, methT),2);
% % for non-speed cells
% bslFR_nonsp = mean(meanFR(nonspeedcell, bslT),2);
% methFR_nonsp = mean(meanFR(nonspeedcell, methT),2);
% % for speed cells
% bslFR_sp = mean(meanFR(speedidx, bslT),2);
% methFR_sp = mean(meanFR(speedidx, methT),2);
% % for negative speed cells
% bslFR_negsp = mean(meanFR(negspeedidx, bslT),2);
% methFR_negsp = mean(meanFR(negspeedidx, methT),2);
% save('FR_groupmean.mat', 'bslFR', 'methFR', 'bslFR_nonsp', 'methFR_nonsp',...
%     'bslFR_sp', 'methFR_sp', 'bslFR_negsp', 'methFR_negsp')
% figure;
% plotSpread([bslFR_nonsp,methFR_nonsp])
% [h,p] = signrank(bslFR_nonsp, methFR_nonsp)
%% plot correlation for nonspeed and speed cells
corr_speed = corrMatrix(:,:,speedidx); corr_speed_blk = corrBlock(:,:,speedidx); 
corr_negspeed = corrMatrix(:,:,negspeedidx); corr_negspeed_blk = corrBlock(:,:,negspeedidx); 
corr_nonspeed = corrMatrix(:,:,nonspeedcell); corr_nonspeed_blk = corrBlock(:,:,nonspeedcell); 
figure;
subplot(2,2,1);
imagesc(nanmean(corrMatrix, 3));
axis square;
title('All Neurons');
subplot(2,2,2);
imagesc(nanmean(corr_nonspeed, 3));
axis square;
title('Non-speed Cells');
subplot(2,2,3);
imagesc(nanmean(corr_speed, 3));
axis square;
title('Speed Cells');
subplot(2,2,4);
imagesc(nanmean(corr_negspeed, 3));
axis square;
title('Negative Speed Cells');
saveas(gcf, 'singleTrialCorr_Speedmod.fig');

figure;
subplot(2,2,1);
imagesc(nanmean(corrBlock, 3));
axis square;
title('All Neurons');
subplot(2,2,2);
imagesc(nanmean(corr_nonspeed_blk, 3));
axis square;
title('Non-speed Cells');
subplot(2,2,3);
imagesc(nanmean(corr_speed_blk, 3));
axis square;
title('Speed Cells');
subplot(2,2,4);
imagesc(nanmean(corr_negspeed_blk, 3));
axis square;
title('Negative Speed Cells');
saveas(gcf, 'trialChunkCorr_Speedmod.fig');

%% plot rate map of non-speed and speed cells(optional)
cellnumber = (1:length(FRS))';
cellnumber = cellnumber(nonspeedcell);
plot_ratemap(FRS, cellID, cellnumber);

%% spatial cell analysis
% identify spatial cells by using spatial information
load(matPath);
load('ratemap_n_corr.mat');
trackLength = 320; %John Wen's track
trial_range = [1:50]; %trials that used to define spatial cells
[spatialcell_sec, spatialcell_spike, avgFR, Tinfo] = ...
    spatial_info_shuffle(S, T, speed, post, posx, sp, trial, trackLength, trial_range);
highfcell = find(avgFR >= 0.1);
% add an additional constraint by using the correlation
half1 = 1:max(trial_range)/2; half2 = max(half1)+1:max(trial_range);
stability = [];
for ii = 1:length(FRS)
    sta = corr2(FRS{ii}(half1,:),FRS{ii}(half2,:));
    stability(ii,1) = sta;
end
stablecell = find(stability > 0.2); %0.2 for MEC, 0.4 for HPC
placecell = intersect(intersect(spatialcell_spike,stablecell),highfcell);
save('placecell.mat', 'placecell','Tinfo','avgFR');
% plot place cells rate map
plot_ratemap(FRS, cellID, placecell');
% calculate firing rate changes of place cells
meanFR_pc = meanFR(placecell,:);
figure;
plot(mean(meanFR_pc))

corrBlock_pc = corrBlock(:,:,placecell);
figure;
imagesc(nanmean(corrBlock_pc,3))
save('placecell.mat', 'meanFR_pc','corrBlock_pc', '-append');

% place cells defined by meth session
trackLength = 320; %John Wen's track
trial_range = methTrials; %trials that used to define spatial cells
[spatialcell_sec, spatialcell_spike, avgFR, Tinfo_meth] = ...
    spatial_info_shuffle(S, T, speed, post, posx, sp, trial, trackLength, trial_range);
save('placecell.mat', 'Tinfo_meth', '-append');

