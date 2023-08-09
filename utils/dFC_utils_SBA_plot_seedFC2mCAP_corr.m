
function dFC_utils_SBA_plot_seedFC2mCAP_corr(file_name, FC_Params)

% seed_name = 'Right SS';
% FC_Params = seed_data.FC_Params;

SpatCorr1M = mean(FC_Params.SpatCorr,1);
SpatCorr1SD = std(FC_Params.SpatCorr,1);

percentile1M = mean(FC_Params.percentile,1);

figure
fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 14;
fig.Position(4) = 9;
fig.Position(1) = 20;
fig.Position(2) = 20;

mseb(percentile1M,SpatCorr1M,SpatCorr1SD);    %This shaded error-bar graphing code is in dFC scripts.

set(gca,'xdir','reverse');
    xlabel('Threshold (percentile)', 'FontSize', 12)
    ylabel('Spatial correlation to FC map (r)', 'FontSize', 12)
    
    
    
    a = get(gca,'XTickLabel');
    set(gca,'XMinorTick','on','XTickLabel',a,'fontsize',12)
%     XMinorTick[on,{off}]
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'fontsize',12)
    
box off
fig.PaperPositionMode = 'auto';


print([file_name], '-dpng','-r600')
    

end