function [posx,post,trial,sp] = reformat_rawdata (savepath, posx, post, trial, sp, datapoint)

% This function aims to keep a subset trials of neuropixel recording data
% with reformat the extended track 640cm into normal track 320 cm
% input args
% posx, post, trial, and sp are all neuropixel raw data
% datapoint: 2 numbers define start and end of the trials to keep, or 4
% numbers define two chunks of data. Currently the code only support 1
% chunk or 2 chunks of data.
% Yanjun Sun, Stanford University, 2/11/2020

% to determine trials to keep
if length(datapoint) == 2
    trials2keep = [datapoint(1):datapoint(2)];
elseif length(datapoint) == 4
    trials2keep = [[datapoint(1):datapoint(2)], [datapoint(3):datapoint(4)]];
end

% change sp
if length(datapoint) == 2
    chunk1 = find(trial == datapoint(1) | trial == datapoint(2));
    t1 = post(chunk1(1));
    t2 = post(chunk1(end));
    spike2keep = sp.st >= t1 & sp.st <=t2;
    sp.st = sp.st(spike2keep);
    sp.clu = sp.clu(spike2keep);
    sp.spikeTemplates = sp.spikeTemplates(spike2keep);
    sp.tempScalingAmps = sp.tempScalingAmps(spike2keep);
elseif length(datapoint) == 4
    chunk1 = find(trial == datapoint(1) | trial == datapoint(2));
    t1 = post(chunk1(1));
    t2 = post(chunk1(end));
    chunk2 = find(trial == datapoint(3) | trial == datapoint(4));
    t3 = post(chunk2(1));
    t4 = post(chunk2(end));
    spike2keep = or(sp.st >= t1 & sp.st <=t2, sp.st >=t3 & sp.st <= t4);
    sp.st = sp.st(spike2keep);
    sp.clu = sp.clu(spike2keep);
    sp.spikeTemplates = sp.spikeTemplates(spike2keep);
    sp.tempScalingAmps = sp.tempScalingAmps(spike2keep);
end

% to determine track length and reformated track length
trackL = max(posx);
trackLr = trackL/2;
overshootf = posx >= trackLr;

% to reformat the extended track 640cm into normal track 320cm
posx(overshootf) = posx(overshootf) - trackLr;

% trials2keep = [[2:51],[101:150]];
% to determine which trials to keep and delete the non-need trial frames
keepf = trial == trials2keep;
keepf = logical(sum(keepf,2));
post(~keepf) = [];
posx(~keepf) = [];

% to make a new vector of trial
t = diff(posx);
td = find(t < -300);
td = [0;td;length(posx)];
trial = NaN(length(posx),1);
n = 1;
for ii = 1:(length(td)-1)
    trial(td(ii)+1:td(ii+1),1) = n;
    n = n+1;
end

% output the data
splitname = strsplit(savepath, '\');
filename = splitname{end};
save(fullfile(savepath, [filename,'.mat']), 'sp', 'posx','post','trial', '-v7.3');

end

