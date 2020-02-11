function [posx_r,post_r,trial_r,sp_r] = reformat_rawdata (posx, post, trial, sp)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


trackL = max(posx);
trackLr = trackL/2;
overshootf = posx >= trackLr;

posx_r = posx;
posx_r(overshootf) = posx(overshootf) - trackLr;

trials2keep = [[31:80],[101:150]];
keepf = trial == trials2keep;
keepf = logical(sum(keepf,2));

post_r = post;
post_r(~keepf) = [];
posx_r(~keepf) = [];



t = diff(posx_r);
td = find(t < -300);
td = [0;td;length(posx_r)];
trial_r = NaN(length(posx_r),1);
n = 1;
for ii = 1:(length(td)-1)
    trial_r(td(ii)+1:td(ii+1),1) = n;
    n = n+1;
end

% change sp
chunk1 = find(trial == 31 | trial == 80);
t1 = post(chunk1(1));
t2 = post(chunk1(end));
chunk2 = find(trial == 101 | trial == 150);
t3 = post(chunk2(1));
t4 = post(chunk2(end));

spike2keep = or(sp.st >= t1 & sp.st <=t2, sp.st >=t3 & sp.st <= t4);
sp_r = sp;
sp_r.st = sp.st(spike2keep);
sp_r.clu = sp.clu(spike2keep);
sp_r.clu = sp.clu(spike2keep);
sp_r.spikeTemplates = sp.spikeTemplates(spike2keep);
sp_r.tempScalingAmps = sp.tempScalingAmps(spike2keep);

save ('F3_190625_MEC_reformat.mat', 'sp', 'posx','post','trial', '-v7.3')









end

