function dFC_utils_plot_imagesc_FC(mat,axC,cnames, cmap)

%mat = matrix to plot using imagesc.
% cmap = the name of the colormap to be used (ex: 'jet')
figure('Position', [15,15,800,750])
imagesc(mat);            % Create a colored plot of the matrix values
colormap(cmap);  % Change the colormap to gray (so higher values are
set(gca, 'XTick', 1:size(mat,1), ...                             % Change the axes tick marks
         'XTickLabel', cnames, ...  %   and tick labels
         'YTick', 1:size(mat,1), ...
         'YTickLabel', cnames, ...
         'TickLength', [0 0]);
 xtickangle(45)
 colorbar
 caxis(axC)     
     
end %functiion