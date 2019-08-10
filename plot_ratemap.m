function plot_ratemap(FRS, cellID, cellnumber, colbar, scale, nperfig, nfigcol)
%% This function is used to plot firing rate map of longitudinally tracked neurons
% Inputs:
% (1) neuronIndividualsf: a cell array of neuron data across different days, each cell is a source2D variable. Obtained by CNMF-E
% (2) behavIndividualsf: a cell array of filtered behavior data, obtained by using Tristan's code
% (3) firinrateAll: a cell array of firing rates of all the neurons across all the recording days
% (4) cellnumber: IDs of the cell that want to be plotted, e.g. 1:size(neuron.trace,1)
% (5) thresh: a vector contains cut off threshold of calcium signals for each neuron, usually defined as 2*SD of neuron.S or neuron.trace
% (6) dataType: 'S' or 'trace'
% (7) scale: max heatmap ploting scale (Hz); use 'auto' if want to plot as the max firing rate of each individual cell
if ~exist('nfigcol', 'var') || isempty(nfigcol);
    nfigcol = 10;
end
if ~exist('nperfig', 'var') || isempty(nperfig);
    nperfig = 30;
end
if ~exist('scale', 'var') || isempty(scale);
    scale = 'auto';
end
if ~exist('colbar', 'var') || isempty(colbar);
    colbar = false;
end
if ~exist('cellnumber', 'var') || isempty(cellnumber);
    cellnumber = 1:length(FRS);
end
%% create a folder to save figures
folderName = 'RatemapFigures';
if ~exist(folderName,'dir')
    mkdir(folderName);
end
fpath=folderName;

Pix_SS = get(0,'screensize');
n = 0;
%% plot rate map figures
for ii = 1:ceil(length(cellnumber)/nperfig)
k = 0;
ax = figure;
set(ax, 'Position', Pix_SS);
hold on;
    % inside each figure
    for jj = (n*nperfig+1) : (n*nperfig+nperfig)
        cellnum = cellnumber(jj);
        k = k+1;
        firingRateSmoothing = FRS{jj};
        x = size(firingRateSmoothing,2);
        y = size(firingRateSmoothing,1);
        maxCount = scale;
        if strcmpi(scale,'auto')
            if isempty(firingRateSmoothing)
                 maxCount = 1;
            else
                 maxCount = max(max(firingRateSmoothing));
            end
        end

        subplot(nperfig/nfigcol,nfigcol,k);
        xlim([1 x]);
        ylim([1 y]);
        hold on
        pcolor(firingRateSmoothing);
        colormap(pink);
        caxis([0,maxCount]);
        shading flat;
        title(['Cell' num2str(cellID(cellnum))]);
        if colbar
            cb = colorbar('eastoutside');
            if strcmpi(scale,'auto')
                set(cb,'Ticks',[0, maxCount]);
            end
        end

        if jj == max(cellnumber)
            break
        end
    end

if strcmpi(scale,'auto');
    saveas(gcf,fullfile(fpath,strcat('Cell',num2str(n*nperfig+1), 'to', num2str(n*nperfig+nperfig),'_ratemap_auto','.fig')));
else
    saveas(gcf,fullfile(fpath,strcat('Cell',num2str(n*nperfig+1), 'to', num2str(n*nperfig+nperfig),'_ratemap','.fig')));
end

n = n+1;
end

end
