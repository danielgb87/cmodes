function dFC_utils_plot_imagesc_matrix_with_values(mat, cmap)

%mat = matrix to plot using imagesc.
% cmap = the name of the colormap to be used (ex: 'jet')
figure('Position', [15,32,800,750])
imagesc(mat);            % Create a colored plot of the matrix values
colormap(cmap);  % Change the colormap to gray (so higher values are
                         %   black and lower values are white)

textStrings = num2str(mat(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:size(mat,1));  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% textColors = repmat(mat(:) > midValue, 1, 3);  % Choose white or black for the
%                                                %   text color of the strings so
%                                                %   they can be easily seen over
%                                                %   the background color
% set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

for c = 1:size(mat,1)
    cnames{c} = ['CAP ' num2str(c)];
end
set(gca, 'XTick', 1:size(mat,1), ...                             % Change the axes tick marks
         'XTickLabel', cnames, ...  %   and tick labels
         'YTick', 1:size(mat,1), ...
         'YTickLabel', cnames, ...
         'TickLength', [0 0]);
 xtickangle(45)
 colorbar
 caxis([-1 1])     
     
end %functiion