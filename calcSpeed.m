function speed = calcSpeed(posx,post)
% function to calculate speed from position and time data
% edited Yanjun Sun 8/11/19
% Malcolm Campbell 5/21/15
% modified 6/6/18 MGC
% all time bins are now equal length
% has not been fully tested yet
%
% inputs:
%     posx: positions
%     post: position time
%     p: params (spatial bin size, etc)
%     smoothingSigma: sigma for smoothing of running speed trace
% outputs:
%     speed: smoothed running speed

%% calculate raw speed
dt = median(diff(post)); 
speed = diff(posx)/dt;

% throw out extreme values and interpolate
speed(speed > 150) = NaN;
speed(speed < -5) = NaN;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = [0;speed];

%% smooth speed trace
% smoothSigma = 10;
% speed = gauss_smoothing(speed,smoothSigma);
speed = smoothts(speed','b',ceil(1/dt));
speed = speed';
end