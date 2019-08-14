function [lick_idx, lickMat, lickAccByTrial] = plot_lick_behav(lickt, lickx, post, posx, trial, binsize)
%plot_lick_behav is to visualize animal's behavior data (lick) while
%running down the linear track in virtual reality
%inputs:
%lickt, lickx, post, trial from neuropixel mat file
%binsize, spatial bin size
%Yanjun Sun, Stanford University, 8/14/2019

if ~exist('binsize', 'var') || isempty(binsize)
    binsize = 5;
end 

folderName = 'LickBehav';
if ~exist(folderName,'dir')
    mkdir(folderName);
end
fpath=folderName;

%% raster plot of lick behavior
[~,~,lick_idx] = histcounts(lickt, post);
figure;
plot(posx(lick_idx),trial(lick_idx),'.', 'color', 'r', 'MarkerSize',7);
xlim([-10, max(posx)+10]);
ylim([0, max(trial)]);
xlabel('Track Length');
ylabel('Trials');
saveas(gcf,fullfile(fpath,'LickRaster.fig'));

%% other plots of lick behavior
licktrial = trial(lick_idx);
edges = 0:binsize:640;
lickMat = [];
for ii = 1:max(trial)
    lickMat(ii,:) = histcounts(lickx(licktrial == ii), edges);
end
% plot licks by spatial bins
figure;
bar(sum(lickMat), 'FaceColor','r','EdgeColor','r');
ylabel('Lick counts');
xlabel('Spatial bins');
saveas(gcf,fullfile(fpath,'LickByPosition.fig'));
% plot lick by trials
figure;
plot(sum(lickMat,2),':');
hold on
plot(sum(lickMat,2),'o', 'color', 'r', 'MarkerSize',4);
xlabel('Trials');
ylabel('Lick counts');
saveas(gcf,fullfile(fpath,'LickByTrial.fig'));

%% plot lick accuracy
% adapted from Kei Masuda
lickAccByTrial = zeros(1,max(trial));
for i = 1:max(trial)
    trialLicks = lickx(trial(lick_idx) == i);
    goodLicks = sum(trialLicks < 10) + sum(trialLicks > max(posx)-10); 
    if trialLicks ~= 0
        lickAccByTrial(i) = goodLicks/numel(trialLicks);
    else
        lickAccByTrial(i) = 0.0;
    end
end
figure;
plot(lickAccByTrial,'o', 'color', 'r', 'MarkerSize',4);
ylim([-0.1, 1.1]);
xlim([-5, max(trial)+5]);
xlabel('Trials');
ylabel('Lick accuracy');
saveas(gcf,fullfile(fpath,'LickAccuracy.fig'));

%% save data
save('lickbehav.mat', 'lick_idx', 'lickMat', 'lickAccByTrial');

end

