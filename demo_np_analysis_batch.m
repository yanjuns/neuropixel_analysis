%% integrate all the results from different mouse
%% define data location
data_loc_mec = {}; % folder path can be paste into cell array
savepath = 'D:\Giocomo_Neuropixels_analysis\NP_results';
save(fullfile(savepath, 'data_location.mat'), 'data_loc_mec', '-append');

%% calculate mean FR
load('data_location.mat')
meanFR_nonsp = []; meanFR_sp = []; meanFR_negsp = [];
for ii = 1:length(data_loc_mec)
    cd(data_loc_mec{ii,1});
    load('FR_groupmean.mat');
    meanFR_nonsp = [meanFR_nonsp; [bslFR_nonsp,methFR_nonsp]];
    meanFR_sp = [meanFR_sp; [bslFR_sp,methFR_sp]];
    meanFR_negsp = [meanFR_negsp; [bslFR_negsp,methFR_negsp]];   
end
save(fullfile(savepath, 'meanFR_groupchange.mat'), 'meanFR_nonsp', 'meanFR_sp', 'meanFR_negsp')

[p1,h1] = signrank(meanFR_nonsp(:,1), meanFR_nonsp(:,2))
[p2,h2] = signrank(meanFR_sp(:,1), meanFR_sp(:,2))
[p3,h3] = signrank(meanFR_negsp(:,1), meanFR_negsp(:,2))
figure;
plotSpread(padcat(meanFR_nonsp(:,1),meanFR_nonsp(:,2),meanFR_sp(:,1),meanFR_sp(:,2),meanFR_negsp(:,1),meanFR_negsp(:,2)),...
    'distributionColors', {[43 57 144]/255, [190 30 45]/255,[43 57 144]/255, [190 30 45]/255,[43 57 144]/255, [190 30 45]/255},...
    'distributionMarkers', '.', 'showMM', 4);
text(1,65, ['P = ', num2str(p1)]);
text(3,65, ['P = ', num2str(p2)]);
text(5,65, ['P = ', num2str(p3)]);
xticklabels({'Baseline','Meth','Baseline','Meth','Baseline','Meth'})
ylabel('Mean Firing Rate (spikes/sec)')
title('Mean firing rate differences before and during MA')
