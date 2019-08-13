%% Neuropixel data analysis streamline
%% define trials for each animal
bslTrials = [1:90]; %specify for each particular mouse 
methTrials = [101:136]; %specify for each particular mouse 
%% Get firing rate map and rate map correlation for each neuron
matPath = 'C:\Users\yanjuns\Desktop\Yanjun Data\Miniscope_CPP_analysis\Neuropixels\E1_190619_johnrepeatingtrack_train100+meth37.mat'
trackLength = 640 %John Wen's extended track
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
ylabel('Speed (cm/s)')

% test if mean firing rate significantly changed
bls = mean(meanFR(:, bslTrials),2);
meth = mean(meanFR(:, methTrials),2);
[p,h] = signrank(bls, meth)
% test speed difference
[h,p] = ttest2(speedbytrial(bslTrials), speedbytrial(methTrials))

%% look at the correlation between firing rate and speed
% overall correlation
corr(speedbytrial(:,bslTrials)', mean(meanFR(:,bslTrials))')
corr(speedbytrial(:,methTrials)', mean(meanFR(:,methTrials))')
% bin by bin correlation
FR_speed_corr = [];
binsize = 10; k = 0;
for ii = 1:binsize:max(trial)
    k = k+1;
    if ii+binsize-1 < max(trial)
        FR_speed_corr(k,:) = corr(speedbytrial(:,ii:ii+9)', mean(meanFR(:,ii:ii+9))');
    else
        FR_speed_corr(k,:) = corr(speedbytrial(:,ii:end)', mean(meanFR(:,ii:end))');
    end
end
[h,p] = ttest2(FR_speed_corr(1:ceil(max(bslTrials)/10),:), FR_speed_corr(11:ceil(max(methTrials)/10),:))


