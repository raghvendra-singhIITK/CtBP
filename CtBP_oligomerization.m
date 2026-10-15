% CtBP1/CtBP2 tetramer assembly as a function of free cytoplasmic NADH
% Cooperative, Hill-equation based model driven by NADH evaluated at six physiologically defined oxygen states.
% Parameters
Kd_CtBP1 = 49.7;
n_CtBP1  = 2.8;
Kd_CtBP2 = 107.1;
n_CtBP2  = 1.6;
free_low  = 0.005;
free      = 0.01;
free_high = 0.015;
stateName  = {'Hyperoxia','Normoxia','Mild hypoxia','Moderate hypoxia', 'Severe hypoxia','Near anoxia'};
NADH_total = [3487; 5224; 10410; 25690; 97330; 219000];
numStates  = numel(stateName);
NADH_free      = free      * NADH_total;
NADH_free_low  = free_low  * NADH_total;
NADH_free_high = free_high * NADH_total;
% Hill equations
hill_function_CtBP1 = @(NADH) (NADH.^n_CtBP1) ./ (Kd_CtBP1^n_CtBP1 + NADH.^n_CtBP1);
hill_function_CtBP2 = @(NADH) (NADH.^n_CtBP2) ./ (Kd_CtBP2^n_CtBP2 + NADH.^n_CtBP2);
Tetramer_Fraction_CtBP1      = hill_function_CtBP1(NADH_free);
Tetramer_Fraction_CtBP1_low  = hill_function_CtBP1(NADH_free_low);
Tetramer_Fraction_CtBP1_high = hill_function_CtBP1(NADH_free_high);
Tetramer_Fraction_CtBP2      = hill_function_CtBP2(NADH_free);
Tetramer_Fraction_CtBP2_low  = hill_function_CtBP2(NADH_free_low);
Tetramer_Fraction_CtBP2_high = hill_function_CtBP2(NADH_free_high);
% Range bars
eb_CtBP1 = Tetramer_Fraction_CtBP1 - Tetramer_Fraction_CtBP1_low;
ea_CtBP1 = Tetramer_Fraction_CtBP1_high - Tetramer_Fraction_CtBP1;
eb_CtBP2 = Tetramer_Fraction_CtBP2 - Tetramer_Fraction_CtBP2_low;
ea_CtBP2 = Tetramer_Fraction_CtBP2_high - Tetramer_Fraction_CtBP2;
fprintf('\nTetramer formed (%%), free NADH = %.1f%% of total cytosolic NADH\n', 100*free);
fprintf('%-18s %10s %10s\n', 'Oxygen state', 'CtBP1', 'CtBP2');
for stateIndex = 1:numStates
    fprintf('%-18s %9.1f%% %9.1f%%\n', stateName{stateIndex}, 100*Tetramer_Fraction_CtBP1(stateIndex), 100*Tetramer_Fraction_CtBP2(stateIndex));
end
fprintf('\n');
% Plot format
axisLineWidth     = 2.5;
tickFontSize      = 26;
axisLabelFontSize = 26;
barWidth          = 1.5;
errorBarLineWidth = 2.5;
errorBarCapSize   = 8;
axesFormat = {'FontName','Arial','FontSize', tickFontSize, 'LineWidth', axisLineWidth, 'XColor', 'k', 'YColor', 'k', 'Box', 'off'};
paralogColor_CtBP1 = [0.1216 0.4667 0.7059];  
paralogColor_CtBP2 = [0.8000 0.4000 0.0000];
pointColor = [0.231 0.490 0.847; 
    0.839 0.153 0.157; 
    0.149 0.651 0.604; 
    0.850 0.650 0.130; 
    0.890 0.467 0.761 
    0.106 0.369 0.125]; 
gridColor   = [0.85 0.85 0.85];
kdLineColor = [0.55 0.55 0.55];
% bar plot
stateNameWrap = {'Hyperoxia','Normoxia','Mild\newlinehypoxia','Moderate\newlinehypoxia','Severe\newlinehypoxia','Near\newlineanoxia'};
barValues     = [Tetramer_Fraction_CtBP1, Tetramer_Fraction_CtBP2];
barErrorBelow = [eb_CtBP1, eb_CtBP2];
barErrorAbove = [ea_CtBP1, ea_CtBP2];
figureBar = figure('Color','w','Units','inches','Position',[1 1 8 4.4]);
axesBar   = axes(figureBar); 
hold(axesBar,'on');
barHandles = bar(axesBar, barValues, 'grouped', 'EdgeColor','none', 'BarWidth', barWidth);
barHandles(1).FaceColor = paralogColor_CtBP1;
barHandles(2).FaceColor = paralogColor_CtBP2;
barXPositions = [barHandles(1).XEndPoints; barHandles(2).XEndPoints]';
errorbar(axesBar, barXPositions, barValues, barErrorBelow, barErrorAbove, 'LineStyle','none', 'Color','k', 'LineWidth', errorBarLineWidth, ...
         'CapSize', errorBarCapSize, 'HandleVisibility','off');
set(axesBar, axesFormat{:}, 'XTick', 1:numStates, 'XTickLabel', stateNameWrap, 'XTickLabelRotation', 0, 'TickLength', [0 0], 'YTick', 0:0.2:1, ...
    'YTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'});
ylim(axesBar, [0 1.1]);
ylabel(axesBar, 'Tetramer fraction', 'FontSize', axisLabelFontSize, 'Color','k');
legend(axesBar, barHandles, {'CtBP1','CtBP2'}, 'Location','northwest', 'Box','off', 'FontSize', tickFontSize);

% hill curve
NADH_curve = logspace(log10(10), log10(5000), 500);
plotHillCurve(NADH_curve, NADH_free, hill_function_CtBP1, Tetramer_Fraction_CtBP1, Kd_CtBP1, paralogColor_CtBP1, pointColor, stateName, ...
              axesFormat, gridColor, kdLineColor, axisLabelFontSize, 'CtBP1');
plotHillCurve(NADH_curve, NADH_free, hill_function_CtBP2, Tetramer_Fraction_CtBP2, Kd_CtBP2, paralogColor_CtBP2, pointColor, stateName, ...
              axesFormat, gridColor, kdLineColor, axisLabelFontSize, 'CtBP2');
function plotHillCurve(NADH_curve, NADH_free, hill_function, Tetramer_Fraction, Kd, curveColor, pointColor, stateName, axesFormat, ...
    gridColor, kdLineColor, axisLabelFontSize, paralogName)
    fig = figure('Color','w','Units','inches','Position',[1 1 6.5 5.3]);
    ax  = axes(fig); 
    hold(ax,'on');
    set(ax, axesFormat{:}, 'XScale','log', 'XLim',[10 5000], 'XTick',[10 100 1000], 'YLim',[0 1.18], 'YTick',0:0.2:1, ...
    'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'}, 'XGrid','on', 'YGrid','on', 'GridColor',gridColor, 'GridAlpha',1, 'GridLineStyle','-');
    plot(ax, [Kd Kd], ax.YLim, '--', 'Color',kdLineColor, 'LineWidth',1);
    text(ax, Kd, 1.10, sprintf('K_d = %.1f nM', Kd), 'Color',kdLineColor, 'FontName','Arial','FontSize',11.5, 'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    plot(ax, NADH_curve, hill_function(NADH_curve), '-', 'Color',curveColor, 'LineWidth',4);
    for i = 1:numel(stateName)
        scatter(ax, NADH_free(i), Tetramer_Fraction(i), 80, pointColor(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.5);
        text(ax, NADH_free(i), Tetramer_Fraction(i), stateName{i}, 'Color','k', 'FontSize',16, 'FontName','Arial', ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
    xlabel(ax, 'Free cytoplasmic NADH (nM)', 'FontSize',axisLabelFontSize, 'Color','k');
    ylabel(ax, 'Tetramer fraction', 'FontSize',axisLabelFontSize, 'Color','k');
    t = title(ax, paralogName, 'FontSize',26, 'FontWeight','normal', 'Color','k');
    t.Units = 'normalized';
    t.Position(1) = 0;
    t.HorizontalAlignment = 'left';
end