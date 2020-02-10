%% integrate all the results from different mouse
%% define data location
data_loc_hpc = {}; % folder path can be paste into cell array
savepath = 'D:\Giocomo_Neuropixels_analysis\NP_results_HPC';
save(fullfile(savepath, 'data_location.mat'), 'data_loc_hpc');

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

%% spatial cell analysis
load('data_location.mat')
meanFRpc = []; corrBlockpc = [];
for ii = 1:length(data_loc_hpc)
    cd(data_loc_hpc{ii,1});
    load('placecell.mat');
    % concatenate meanFR_pc
    meanFRpc = padconcatenation(meanFRpc,meanFR_pc,1);
    % concatenate corrBlock_pc
    if size(corrBlock_pc,1) < 20
       temp =  NaN(20,20,size(corrBlock_pc,3));
       temp(1:size(corrBlock_pc,1),1:size(corrBlock_pc,2),1:size(corrBlock_pc,3)) = ...
           corrBlock_pc;
       corrBlock_pc = temp;
    end
    corrBlockpc = cat(3, corrBlockpc, corrBlock_pc);
end
meanFRall = []; speedall = [];
for ii = 1:length(data_loc_hpc)
    cd(data_loc_hpc{ii,1});
    load('ratemap_n_corr.mat');
    % concatenate meanFR_pc
    meanFRall = padconcatenation(meanFRall,meanFR,1);
    speedall = padconcatenation(speedall,speedbytrial,1);
end

save(fullfile(savepath, 'spatialcell_analysis.mat'), 'meanFRpc','corrBlockpc','meanFRall','speedall')

% plot the results of firing rate and speed
frbslall = nanmean(meanFRall(:,1:100),2);
frmethall = nanmean(meanFRall(:,101:200),2);
[p1,h1] = signrank(frbslall,frmethall)
frbsl = nanmean(meanFRpc(:,1:100),2);
frmeth = nanmean(meanFRpc(:,101:200),2);
[p2,h2] = signrank(frbsl,frmeth)
frbslsp = nanmean(speedall(:,1:100),2);
frmethsp = nanmean(speedall(:,101:200),2);
[p3,h3] = signrank(frbslsp,frmethsp)

figure;
subplot(3,1,1)
stdshade(meanFRall, 0.25, 'b');
ylim([15 30]);
ylabel('Avg firing rate (spike/sec)');
title('Firing rate change across all neurons');
text(140,26, ['P = ', num2str(p1)]);
subplot(3,1,2)
stdshade(meanFRpc, 0.25);
xlabel('Trial');
ylabel('Avg firing rate (spike/sec)');
title('Firing rate change of spatial cells');
text(140,40, ['P = ', num2str(p2)]);
subplot(3,1,3)
stdshade(speedall, 0.25, 'g');
xlabel('Trial');
ylabel('cm / sec');
title('Change of speed');
text(140,35, ['P = ', num2str(p3)]);

% plot the results of correlation
figure;
imagesc(nanmean(corrBlockpc,3));
axis image
colorbar
xlabel('Trial block (10 trials avg / block)')
ylabel('Trial block (10 trials avg / block)')

%speed effect
speedcell = cellID(speedidx);
intersect(placecell, speedcell)
%there is ovreall very rare speed modulation of speed cells

