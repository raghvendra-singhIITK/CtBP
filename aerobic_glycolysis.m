% Kinetic model for aerobic glycolysis
function aerobic_glycolysis
% Different oxygen fraction levels for simulation
oxy_levels = [0.02 0.05 0.2 0.5 1 1.5]; % near-anoxia to hyperoxia
tspan = [0 10000];
num_oxy = numel(oxy_levels);
[species, y0, species_count] = Model_species();
parameter = Model_parameters();
% figure formatting
set(0, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontSize', 26, 'DefaultTextFontName', 'Times New Roman', ...
   'DefaultTextFontSize', 26, 'DefaultLineLineWidth', 3, 'DefaultAxesLineWidth', 1.5, 'DefaultAxesTickLength', [0.01 0.01], ...
   'DefaultFigureColor', 'w')
idx_NADH = 17;
% Simulations
all_t = cell(1, num_oxy);
all_conc = cell(1, num_oxy);
opts = odeset('NonNegative', 1:species_count);
for k = 1:num_oxy
    [all_t{k}, all_conc{k}] = ode15s(@(t, y) Model(t, y, parameter, oxy_levels(k)), tspan, y0, opts);
    NADH_ss = all_conc{k}(end, 17);
    fprintf('%g oxy = NADH %.4g mM\n', oxy_levels(k), NADH_ss);
end

% NADH bar plot
free_fraction = 0.01;
col  = [149 188 241; 38 83 146] / 255; % bar colours
[~, ord] = sort(oxy_levels, 'descend');
tot  = 1000 * cellfun(@(c) c(end,idx_NADH), all_conc(ord));
tot(tot <= 0) = NaN;
vals = [tot; free_fraction*tot];
y_log = 10^(log10(min(vals(2,:))) - 0.5);
y_hieght = 10^(log10(max(vals(1,:))) + 0.4);
xlab = {'Hyperoxia','Normoxia','Mild\newlinehypoxia','Moderate\newlinehypoxia','Severe\newlinehypoxia','Near\newlineanoxia'};

if numel(xlab) ~= num_oxy, xlab = arrayfun(@(x) sprintf('%.2f', x), oxy_levels(ord), 'UniformOutput', false); 
end
fig = figure('Units','inches', 'Position',[1 1 8 6.46], 'PaperUnits','inches', 'PaperSize',[8 6.46], ...
             'PaperPosition',[0 0 8 6.46], 'Color','w');
ax  = axes('Parent',fig, 'NextPlot','add');
for r = 1:2
    b(r) = bar(ax, 1:num_oxy, vals(r,:), 0.55, 'FaceColor',col(r,:), 'EdgeColor','none', 'BaseValue',y_log);
end
set(ax, 'YScale','log', 'YLim',[y_log y_hieght], 'XLim',[0.45 num_oxy+0.55], 'XTick',1:num_oxy, 'XTickLabel',xlab, ...
        'YMinorTick','on', 'TickDir','out', 'TickLength',[0.0125 0.0125], 'Box','off', 'Layer','top', ...
        'LineWidth', 2.5, 'FontName','Arial', 'FontSize',26, 'Position',[0.152 0.129 0.842 0.818])
ax.XAxis.TickLength = [0 0];
ylabel(ax, ['Cytosolic NADH (' char(181) 'M)'], 'FontName','Arial', 'FontSize',30)
lg = legend(ax, b, {'Total NADH', sprintf('Free NADH (%g%%)', 100*free_fraction)}, ...
            'Location','northwest', 'Box','off', 'FontName','Arial', 'FontSize',26);
           
end

% Species involved in the model
% These 21 species represent the glycolytic intermediates
function [species, y0, species_count] = Model_species()
species = {'Gluout', 'Gluin', 'G6P', 'F6P', 'F16BP', 'DHAP', 'G3P', ...
    'BPG', 'PG3', 'PG2', 'PEP', 'PYR', 'LACin', 'ATP', 'ADP', 'NAD', ... 
    'NADH', 'Pi', 'GLY', 'AMP', 'LACout'}; 
% Initial concentrations of these species in mM
% Values as defined in Marín-Hernández et al. 2011
Gluout = 5;
Gluin = 6.2;
G6P = 5.3;
F6P = 1.5;
F16BP = 25;
DHAP = 10;
G3P = 0.9;
BPG = 0.01; 
PG3 = 0.01;
PG2 = 0.01;
PEP = 0.1;
PYR = 2.1;
LACin = 27;
ATP = 5.6;
ADP = 2.4;
NAD = 1.3;
NADH = 0.01;
Pi = 4.8;
GLY = 26;
AMP = 3.3;
LACout = 0.10;

y0 = [Gluout, Gluin, G6P, F6P, F16BP, DHAP, G3P, BPG, PG3, PG2, PEP, PYR, LACin, ATP, ADP, NAD, NADH, Pi, GLY, AMP, LACout];
species_count = length(species);
end

% Kinetic parameters involved in the model
% Reference for the parameter values used are given in Table S2
function parameter = Model_parameters()
parameter = struct();
parameter.GLUT_Vmf = 0.055;
parameter.GLUT_Keq = 1;
parameter.GLUT_KmGluout = 0.52;
parameter.GLUT_KmGluin = 0.52;
parameter.HK_Vmf = 0.46;
parameter.HK_Keq = 651;
parameter.HK_KmGlu = 0.18;
parameter.HK_KmATP = 0.99;
parameter.HK_KmADP = 3.5;
parameter.HK_KiG6P = 0.02;
parameter.HPI_Vmf = 4.9;
parameter.HPI_Vmr = 3.4;
parameter.HPI_KmG6P = 0.9;
parameter.HPI_KmF6P = 0.07;
parameter.HPI_KiERY4P = 0.0017;
parameter.HPI_KiF16BP = 0.17;
parameter.HPI_Ki6PG = 0.0094;
parameter.PFK_Vmf = 0.273;
parameter.PFK_Keq = 247;
parameter.PFK_KmF6P = 4.6;
parameter.PFK_KmATP = 0.048;
parameter.PFK_KmF16BP = 5;
parameter.PFK_KmADP = 5;
parameter.PFK_KiATP = 1.75;
parameter.PFK_KiCIT = 3.9;
parameter.PFK_KaF26BP = 0.00018;
parameter.PFK_BETA = 2.35;
parameter.PFK_ALPHA = 4.47;
parameter.PFK_L = 13;
parameter.ALDO_Vmf = 0.23;
parameter.ALDO_Vmr = 0.18;
parameter.ALDO_KmF16BP = 0.01;
parameter.ALDO_KmG3P = 0.16;
parameter.ALDO_KmDHAP = 0.08;
parameter.TPI_Vmf = 5.6;
parameter.TPI_Vmr = 56;
parameter.TPI_KmDHAP = 1.9;
parameter.TPI_KmG3P = 0.41;
parameter.GAPDH_Vmf = 1;
parameter.GAPDH_Vmr = 0.9;
parameter.GAPDH_KmG3P = 0.29;
parameter.GAPDH_KmBPG = 0.02;
parameter.GAPDH_KmNAD = 0.08;
parameter.GAPDH_KmNADH = 0.004;
parameter.GAPDH_KmPi = 11;
parameter.PGK_Vmf = 27;
parameter.PGK_Vmr = 4.3;
parameter.PGK_KmBPG = 0.035;
parameter.PGK_Km3PG = 0.12;
parameter.PGK_KmADP = 0.67;
parameter.PGK_KmATP = 0.15;
parameter.PGK_ALPHA = 1;
parameter.PGK_BETA = 1;
parameter.PGAM_Vmf = 20;
parameter.PGAM_Vmr = 1.3;
parameter.PGAM_Km3PG = 0.18;
parameter.PGAM_Km2PG = 0.04;
parameter.ENO_Vmf = 0.51;
parameter.ENO_Vmr = 0.74;
parameter.ENO_Km2PG = 0.16;
parameter.ENO_KmPEP = 0.04;
parameter.PYK_Vmf = 6.6;
parameter.PYK_Keq = 195172;
parameter.PYK_KmPEP = 0.4;
parameter.PYK_KmADP = 0.3;
parameter.PYK_KmPYR = 10;
parameter.PYK_KmATP = 0.86;
parameter.PYK_KaF16BP = 0.0004;
parameter.PYK_KiATP = 2.5;
parameter.PYK_L = 1;
parameter.LDH_Vmf = 13.4;
parameter.LDH_Vmr = 1.8;
parameter.LDH_KmPYR = 0.13;
parameter.LDH_KmLAC = 4.7;
parameter.LDH_KmNAD = 0.07;
parameter.LDH_KmNADH = 0.002;
parameter.LDH_ALPHA = 1;
parameter.LDH_BETA = 1;
parameter.MCT_Vmf = 0.68;
parameter.MCT_Keq = 1.99;
parameter.MCT_KmLAC = 4.68;
parameter.AS_K = 0.0042;
parameter.AK_Kf = 1;
parameter.AK_Kb = 2.26;
parameter.DH_Kf = 250;
parameter.DH_Kr = 1;
parameter.MPM_Vmf = 1.8;
parameter.MPM_KmPYR = 0.9;
parameter.MPM_KmADP = 0.022;
parameter.MPM_KmPi = 4.0;
parameter.MPM_Kmoxy = 0.0025;
end 

% Main model function for the aerobic glycolysis
% Defines the 21 ODEs, where the reaction velocities are defined based on
% each enzymes' catalytic mechanism. also stated in Table S1
function dydt = Model(t, y, parameter, oxy)
% Loading species and parameters
[species, ~, ~] = Model_species();
for i=1:length(species)
    eval(sprintf('%s = y(%d);', species{i}, i));
end
parameter_names = fieldnames(parameter);
for j=1:length(parameter_names)
    eval(sprintf('%s = parameter.%s;', parameter_names{j}, parameter_names{j}));
end
% Regulatory species with fixed concentration
ERY4P = 1;
PG6 = 0.35;
F26BP = 0.006;
CIT = 1.7;
% Velocity equation for each enzyme
v_GLUT = @(Gluout, Gluin) ((GLUT_Vmf*(Gluout-(Gluin/GLUT_Keq)))/(GLUT_KmGluout*(1+(Gluin/GLUT_KmGluin))+Gluout));
v_HK = @(Gluin, ATP, G6P, ADP) (((HK_Vmf/(HK_KmGlu*HK_KmATP))*((Gluin*ATP)-(G6P*ADP/HK_Keq)))/(1+(Gluin/HK_KmGlu)+(ATP/HK_KmATP)+(Gluin*ATP/(HK_KmGlu*HK_KmATP))+(G6P/HK_KiG6P)+(ADP/HK_KmADP)+(G6P*ADP/(HK_KiG6P*HK_KmADP))+(Gluin*ADP/(HK_KmGlu*HK_KmADP))+(ATP*G6P/(HK_KmATP*HK_KiG6P))));
v_HPI = @(G6P, F6P, ERY4P, PG6, F16BP) (((HPI_Vmf*G6P/HPI_KmG6P)-(HPI_Vmr*F6P/HPI_KmF6P))/(1+(G6P/HPI_KmG6P)+(F6P/HPI_KmF6P)+(ERY4P/HPI_KiERY4P)+(PG6/HPI_Ki6PG)+(F16BP/HPI_KiF16BP)));
v_PFK = @(ATP, F26BP, F6P, CIT, ADP, F16BP) (PFK_Vmf* ( ( (ATP/PFK_KmATP)/(1 + ATP/PFK_KmATP))*((1 + (PFK_BETA*F26BP)/(PFK_ALPHA*PFK_KaF26BP))/(1+(F26BP)/(PFK_ALPHA*PFK_KaF26BP)) )*(( (F6P*(1 + F26BP/(PFK_ALPHA*PFK_KaF26BP)))/(PFK_KmF6P*(1 + F26BP/PFK_KaF26BP)))*  ( 1+(F6P*(1 + F26BP/(PFK_ALPHA*PFK_KaF26BP))) /(PFK_KmF6P*(1 + F26BP/PFK_KaF26BP)) )^3 ) /( (PFK_L*(1 + CIT/PFK_KiCIT)^4*(1 + ATP/PFK_KiATP)^4) /(1 + F26BP/PFK_KaF26BP)^4 + ( 1 + (F6P*(1 + F26BP/(PFK_ALPHA*PFK_KaF26BP))) / (PFK_KmF6P*(1 + F26BP/PFK_KaF26BP)) )^4 )  - ( (ADP*F16BP)/(PFK_KmADP*PFK_KmF16BP*PFK_Keq) ) /( 1 + ADP/PFK_KmADP + F16BP/PFK_KmF16BP +  (ADP*F16BP)/(PFK_KmADP*PFK_KmF16BP) )  ));
v_ALDO = @(F16BP, DHAP, G3P) (((ALDO_Vmf*F16BP/ALDO_KmF16BP)-(ALDO_Vmr*(DHAP/ALDO_KmDHAP)*(G3P/ALDO_KmG3P)))/(1+(F16BP/ALDO_KmF16BP)+(DHAP/ALDO_KmDHAP)+(G3P/ALDO_KmG3P)+((G3P/ALDO_KmG3P)*(DHAP/ALDO_KmDHAP))));
v_TPI = @(DHAP, G3P) (((TPI_Vmf*DHAP/TPI_KmDHAP)-(TPI_Vmr*G3P/TPI_KmG3P))/(1+(DHAP/TPI_KmDHAP)+(G3P/TPI_KmG3P)));
v_GAPDH = @(NAD, G3P, Pi, BPG, NADH) (((GAPDH_Vmf*NAD*G3P*Pi/(GAPDH_KmNAD*GAPDH_KmG3P*GAPDH_KmPi))-(GAPDH_Vmr*BPG*NADH/(GAPDH_KmNADH*GAPDH_KmBPG)))/(1+(NAD/GAPDH_KmNAD)+((NAD/GAPDH_KmNAD)*(G3P/GAPDH_KmG3P))+((NAD/GAPDH_KmNAD)*(G3P/GAPDH_KmG3P)*(Pi/GAPDH_KmPi))+((BPG/GAPDH_KmBPG)*(NADH/GAPDH_KmNADH))+(NADH/GAPDH_KmNADH)));
v_PGK = @(BPG, ADP, PG3, ATP) ((((PGK_Vmf*BPG*ADP)/(PGK_ALPHA*PGK_KmBPG*PGK_KmADP))-((PGK_Vmr*PG3*ATP)/(PGK_BETA*PGK_Km3PG*PGK_KmATP)))/(1+(BPG/PGK_KmBPG)+(ADP/PGK_KmADP)+((BPG*ADP)/(PGK_ALPHA*PGK_KmBPG*PGK_KmADP))+(PG3/PGK_Km3PG)+(ATP/PGK_KmATP)+((PG3*ATP)/(PGK_BETA*PGK_Km3PG*PGK_KmATP))));
v_PGAM = @(PG3, PG2) (((PGAM_Vmf*PG3/PGAM_Km3PG)-(PGAM_Vmr*PG2/PGAM_Km2PG))/(1+(PG3/PGAM_Km3PG)+(PG2/PGAM_Km2PG)));
v_ENO = @(PG2, PEP) (((ENO_Vmf*PG2/ENO_Km2PG)-(ENO_Vmr*PEP/ENO_KmPEP))/(1+(PEP/ENO_KmPEP)+(PG2/ENO_Km2PG)));
v_PYK = @(ADP, PEP, ATP, F16BP, PYR) (PYK_Vmf*(( (ADP/PYK_KmADP) / (1 + ADP/PYK_KmADP))*( (PEP/PYK_KmPEP) * (1 + PEP/PYK_KmPEP)^3 ) /( (PYK_L*(1 + ATP/PYK_KiATP)^4) /  (1 + F16BP/PYK_KaF16BP)^4 +   (1 + PEP/PYK_KmPEP)^4 )    -  ( (ATP*PYR)/(PYK_KmATP*PYK_KmPYR*PYK_Keq))/(1+ATP/PYK_KmATP + PYR/PYK_KmPYR + (ATP*PYR)/(PYK_KmATP*PYK_KmPYR) )));
v_LDH = @(NADH, PYR, LACin, NAD) ((((LDH_Vmf*NADH*PYR)/(LDH_ALPHA*LDH_KmNADH*LDH_KmPYR))-((LDH_Vmr*LACin*NAD)/(LDH_BETA*LDH_KmLAC*LDH_KmNAD)))/(1+(NADH/LDH_KmNADH)+(PYR/LDH_KmPYR)+((NADH*PYR)/(LDH_ALPHA*LDH_KmNADH*LDH_KmPYR))+(LACin/LDH_KmLAC)+(NAD/LDH_KmNAD)+((LACin*NAD)/(LDH_BETA*LDH_KmLAC*LDH_KmNAD))));
v_AS = @(ATP) (AS_K*ATP);
v_AK = @(ATP,AMP,ADP) (AK_Kf*ATP*AMP) - (AK_Kb*(ADP^2));
v_DH = @(NADH, NAD, oxy) (DH_Kf * NADH* oxy) -(DH_Kr*NAD);
v_MCT = @(LACin, LACout) ((MCT_Vmf*(LACin-(LACout/MCT_Keq)))/(MCT_KmLAC+LACin+LACout));
v_MPM = @(PYR, ADP, Pi, oxy) MPM_Vmf*(PYR/(MPM_KmPYR + PYR))* (ADP/(MPM_KmADP + ADP))* (Pi/ (MPM_KmPi  + Pi ))* (oxy/ (MPM_Kmoxy + oxy));
% Fixed fluxes for the interacting pathway branches
v_GD = 1.2*10^-3;
v_GS = 2.2*10^-3;
v_PPPB = 9.5*10^-5;
% ODES
dydt = zeros(21, 1);
%gluout
dydt(1) = -v_GLUT(Gluout, Gluin); 
%gluin
dydt(2) = v_GLUT(Gluout, Gluin) - v_HK(Gluin, ATP, G6P, ADP);
%g6p
dydt(3) = v_HK(Gluin, ATP, G6P, ADP) - v_HPI(G6P, F6P, ERY4P, PG6, F16BP) - v_PPPB - v_GS + v_GD;
%f6p
dydt(4) = v_HPI(G6P, F6P, ERY4P, PG6, F16BP) - v_PFK(ATP, F26BP, F6P, CIT, ADP, F16BP) + (v_PPPB*(2/3));
%f16p
dydt(5) = v_PFK(ATP, F26BP, F6P, CIT, ADP, F16BP) - v_ALDO(F16BP, DHAP, G3P);
%dhap
dydt(6) = v_ALDO(F16BP, DHAP, G3P) - v_TPI(DHAP, G3P);
%g3p
dydt(7) = v_TPI(DHAP, G3P) + v_ALDO(F16BP, DHAP, G3P) - v_GAPDH(NAD, G3P, Pi, BPG, NADH) + (v_PPPB*(1/3));
%bpg
dydt(8) = v_GAPDH(NAD, G3P, Pi, BPG, NADH) - v_PGK(BPG, ADP, PG3, ATP);
%3pg
dydt(9) = v_PGK(BPG, ADP, PG3, ATP) - v_PGAM(PG3, PG2);
%2pg
dydt(10) = v_PGAM(PG3, PG2) - v_ENO(PG2, PEP);
%pep
dydt(11) = v_ENO(PG2, PEP) - v_PYK(ADP, PEP, ATP, F16BP, PYR);
%pyr
dydt(12) = v_PYK(ADP, PEP, ATP, F16BP, PYR) - v_LDH(NADH, PYR, LACin, NAD) - v_MPM(PYR, ADP, Pi, oxy);
%lacin
dydt(13) = v_LDH(NADH, PYR, LACin, NAD) - v_MCT(LACin, LACout);
%atp
dydt(14) = -v_HK(Gluin, ATP, G6P, ADP) - v_PFK(ATP, F26BP, F6P, CIT, ADP, F16BP) + v_PGK(BPG, ADP, PG3, ATP) + v_PYK(ADP, PEP, ATP, F16BP, PYR) - v_GS - v_AS(ATP) - v_AK(ATP,AMP,ADP) + (12.5*v_MPM(PYR, ADP, Pi, oxy));
%adp
dydt(15) = v_HK(Gluin, ATP, G6P, ADP) + v_PFK(ATP, F26BP, F6P, CIT, ADP, F16BP) - v_PGK(BPG, ADP, PG3, ATP) - v_PYK(ADP, PEP, ATP, F16BP, PYR) + v_GS + v_AS(ATP) + (2*v_AK(ATP,AMP,ADP)) - (12.5*v_MPM(PYR, ADP, Pi, oxy));
%nad
dydt(16) = -v_GAPDH(NAD, G3P, Pi, BPG, NADH) + v_LDH(NADH, PYR, LACin, NAD) + v_DH(NADH,NAD, oxy) ; 
%nadh
dydt(17) = v_GAPDH(NAD, G3P, Pi, BPG, NADH) - v_LDH(NADH, PYR, LACin, NAD) - v_DH(NADH, NAD, oxy) ; 
%pi
dydt(18) = -v_GAPDH(NAD, G3P, Pi, BPG, NADH) + v_AS(ATP) + (2*v_GS) - v_GD - (12.5*v_MPM(PYR, ADP, Pi, oxy));
%gly
dydt(19) = v_GS - v_GD;
%amp
dydt(20) = -v_AK(ATP,AMP,ADP);
% lacout
dydt(21) = v_MCT(LACin, LACout);
end