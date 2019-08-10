%% Neuropixel data analysis streamline
%% Get firing rate map and rate map correlation for each neuron
matPath = 'C:\Users\yanjuns\Desktop\Yanjun Data\Miniscope_CPP_analysis\Neuropixels\E1_190619_johnrepeatingtrack_train100+meth37.mat'
trackLength = 640 %John Wen's extended track
trials_per_block = 10;
[S, T, FR, FRS, corrMatrix, S_block, T_block, FR_block, FRS_block, corrBlock, cellID]...
    = get_ratemap_n_corr(matPath, trackLength, trials_per_block);
% get mean firing rate for each cell and each trial
meanFR = NaN(size(S,2), size(T,1));
for ii = 1:size(S,2);
    for jj = 1:size(T,1);
        meanFR(ii,jj) = sum(S{ii}(jj,:))/sum(T(jj,:));
    end
end
save('ratemap_n_corr.mat', 'S', 'T', 'FR', 'FRS', 'corrMatrix',...
    'S_block', 'T_block', 'FR_block', 'FRS_block', 'corrBlock', 'cellID','meanFR', '-v7.3');
%% plot ratemap figures
plot_ratemap(FRS, cellID)
%print figure into eps if needed
figname = 'Cell1to30_ratemap_auto';
fig = openfig(figname);
print (fig, '-painters', '-depsc', [figname, '.png']);
%plot a specific one if needed
figure;
cellnum = find(cellID == 170);
pcolor(FRS{cellnum});
shading flat
colormap pink

%% plot ratemap and mean firing rate
figure;
plot(std(meanFR));
xlabel('Trials');
ylabel('Avg Firing Rate of all neurons (Hz)');

% test if mean firing rate significantly changed
blsTrials = [1:90]; %specify for each particular mouse 
methTrials = [101:136]; %specify for each particular mouse 
bls = mean(meanFR(:, blsTrials),2);
meth = mean(meanFR(:, methTrials),2);
[p,h] = signrank(bls, meth)
