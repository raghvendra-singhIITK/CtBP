% CtBP1/CtBP2 utilization hypotheses
% Parameters
Kd_CtBP1 = 49.7;
n_CtBP1 = 2.8;
Kd_CtBP2 = 107.1;
n_CtBP2 = 1.6;
N_HIGH = 2;
stateName  = {'Hyperoxia','Normoxia','Mild hypoxia','Moderate hypoxia','Severe hypoxia','Near anoxia'};
NADH_total = [3487; 5224; 10410; 25690; 97330; 219000];
NADH_free  = 0.01 * NADH_total;
idxNormoxia = find(strcmp(stateName,'Normoxia'));
idxModHypoxia   = find(strcmp(stateName,'Moderate hypoxia'));
idxSevHyppoxia   = find(strcmp(stateName,'Severe hypoxia'));
idxNearAnoxia = find(strcmp(stateName,'Near anoxia'));
x_plot = linspace(0.01, 1100, 4000);
stateColors = [0.231 0.490 0.847;
               0.839 0.153 0.157;
               0.149 0.651 0.604;
               0.850 0.650 0.130;
               0.890 0.467 0.761;
               0.106 0.369 0.125];
paralogColor_CtBP1 = [0.1216 0.4667 0.7059];
paralogColor_CtBP2 = [0.8000 0.4000 0.0000];
grayLinear = [0.65 0.65 0.65];
grayMM     = [0.45 0.45 0.45];
grayHill   = [0.25 0.25 0.25];
offDx_CtBP1 = [  0,  18, -18,  70, -340, 0];
offDy_CtBP1 = [-11, -15,  10,  16,  -12, 0];
offHa_CtBP1 = {'center','left','right','left','right','center'};
offDx_CtBP2 = [  0, -18,  55, -50, -290, 0];
offDy_CtBP2 = [ 10, -14,   9,  13,  -14, 0];
offHa_CtBP2 = {'center','right','left','right','right','center'};
axisLineWidth     = 2.5;
tickFontSize      = 26;
axisLabelFontSize = 26;
axesFormat = {'FontName','Arial','FontSize',tickFontSize,'LineWidth',axisLineWidth,'XColor','k','YColor','k','Box','on'};
plotUtilizationCurve(NADH_free, Kd_CtBP1, n_CtBP1, N_HIGH, x_plot, stateName, stateColors, paralogColor_CtBP1, grayLinear, grayMM, grayHill, ...
    offDx_CtBP1, offDy_CtBP1, offHa_CtBP1, axesFormat, axisLabelFontSize, idxNormoxia, idxModHypoxia, idxSevHyppoxia, idxNearAnoxia, 'CtBP1');
plotUtilizationCurve(NADH_free, Kd_CtBP2, n_CtBP2, N_HIGH, x_plot, stateName, stateColors, paralogColor_CtBP2, grayLinear, grayMM, grayHill, ...
    offDx_CtBP2, offDy_CtBP2, offHa_CtBP2, axesFormat, axisLabelFontSize, idxNormoxia, idxModHypoxia, idxSevHyppoxia, idxNearAnoxia, 'CtBP2');
function plotUtilizationCurve(NADH_free, Kd, n, N_HIGH, x_plot, stateName, stateColors, curveColor, grayLinear, grayMM, grayCoop, ...
    offDx, offDy, offHa, axesFormat, axisLabelFontSize, idxNormoxia, idxModHyp, idxSevHyp, idxNearAnox, paralogName)
    x1 = NADH_free(idxNormoxia);
    x2 = NADH_free(idxModHyp);
    y1 = Hill_Fn(x1, Kd, n);
    y2 = Hill_Fn(x2, Kd, n);
    Line_Fn = @(x) y1 + (y2 - y1)/(x2 - x1) * (x - x1);
    [Vmax_mm, Km_mm] = solveFit(x1, y1, x2, y2, 1);
    [Vmax_hill, K_hill]  = solveFit(x1, y1, x2, y2, N_HIGH);
    fig = figure('Color','w','Units','inches','Position',[1 1 7.6 6.2]);
    ax  = axes(fig); hold(ax,'on');
    set(ax, axesFormat{:}, 'Box','off');
    yline(ax, 100, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',1.5);
    text(ax, 15, 104, '100%', 'Color',[0.55 0.55 0.55], 'FontName','Arial', 'FontSize',16, 'HorizontalAlignment','left');
    hLinear = plot(ax, x_plot, Line_Fn(x_plot), '--', 'Color',grayLinear, 'LineWidth',2.3);
    hMM = plot(ax, x_plot, MM_Fn(x_plot, Vmax_mm, Km_mm), ':', 'Color', grayMM, 'LineWidth',2.3);
    hHill = plot(ax, x_plot, Hill_Fn(x_plot, K_hill, N_HIGH, Vmax_hill),'-.', 'Color', grayCoop, 'LineWidth',2.3);
    hActual = plot(ax, x_plot, Hill_Fn(x_plot, Kd, n), '-', 'Color',curveColor, 'LineWidth',4);
    for i = 1:numel(stateName)
        xf = NADH_free(i);
        if xf <= 1100
            yv = Hill_Fn(xf, Kd, n);
            scatter(ax, xf, yv, 80, stateColors(i,:), 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.5);
            text(ax, xf + offDx(i), yv + offDy(i), stateName{i}, 'Color','k', 'FontName','Arial', 'FontSize',16, 'HorizontalAlignment',offHa{i});
        end
    end
    xlim(ax, [0 1100]); ylim(ax, [0 140]);
    xlabel(ax, 'Free cytoplasmic NADH (nM)', 'FontName','Arial', 'FontSize',axisLabelFontSize, 'Color','k');
    ylabel(ax, sprintf('%s fraction (%%)', paralogName), 'FontName','Arial', 'FontSize',axisLabelFontSize, 'Color','k');
    t = title(ax, paralogName, 'FontSize',26, 'FontWeight','normal', 'Color','k');
    t.Units = 'normalized';
    t.Position(1) = 0;
    t.HorizontalAlignment = 'left';
    legend(ax, [hLinear hMM hHill hActual], {'Linear utilization', 'Michaelis–Menten utilization (n = 1)', 'Cooperative utilization (n = 2)', ...
        sprintf('Production (Hill, n = %.1f, K_d = %.1f nM)', n, Kd)}, 'Location','southeast', 'FontName','Arial', 'FontSize',14, 'Box','off', 'Color','k');
    xsh = NADH_free(idxSevHyp);
    xna = NADH_free(idxNearAnox);
    fprintf('%s off-chart extrapolation:\n', paralogName);
    fprintf('Severe Hypoxia (%.0f nM): line = %.0f%%, MM = %.0f%%\n', xsh, Line_Fn(xsh), MM_Fn(xsh, Vmax_mm, Km_mm));
    fprintf('Near Anoxia (%.0f nM): line = %.0f%%, MM = %.0f%%\n\n', xna, Line_Fn(xna), MM_Fn(xna, Vmax_mm, Km_mm));
end
function y = Hill_Fn(x, K, n, Vmax)
    if nargin < 4, Vmax = 100; end
    y = Vmax .* x.^n ./ (K.^n + x.^n);
end
function y = MM_Fn(x, Vmax, Km)
    y = Vmax .* x ./ (Km + x);
end
function [Vmax, K] = solveFit(x1, y1, x2, y2, n)
    a = x1^n; b = x2^n;
    Vmax = y1*y2*(a - b) / (a*y2 - b*y1);
    K = (a*(Vmax - y1)/y1)^(1/n);
end