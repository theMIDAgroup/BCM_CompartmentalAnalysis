clc
clear
close all

set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

%% Define colors
dark_red = [0.8 0 0];
dark_green = [0 0.8 0];

%% Define paths and load data
path_results = './results/';
path_figures = './figures/';

load(fullfile(path_results, 'K_#0 CT26 NO STS 2012.11.13 PRIMA PET.mat'));
Ca_minus = Ca(t) - 2*sqrt(Ca(t)); err_Ca = abs(Ca(t) - Ca_minus);
Ct_minus = Ct - sqrt(Ct); err_Ct = abs(Ct - Ct_minus);

%% Time-dependent ROI concentrationcurve of the CT26 tumor C_T and its 
%% standard deviation, related to experiment m1. 
figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf, 'Color','none');
errorbar(t,Ct,err_Ct,'Color',dark_green,...
    'LineStyle','-','Marker','o','LineWidth',3,'Markersize',5);
axis([0 40 0 350]); 
% axis square; grid on;
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',32,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',32,'Interpreter','Latex');
% title({'$\mathcal{C}_T$'},'FontSize',25,'Interpreter','Latex')

fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'CT_std'),'-dpng','-r300')

%% Time-dependent concentration curve of the arterial input function C_i 
%% and itsstandard deviation, related to experiment m1
figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf, 'Color','none');
errorbar(t,Ca(t),err_Ca,'Color',dark_red,...
    'LineStyle','-','Marker','o','LineWidth',3,'Markersize',5);
axis([0 40 0 1600]); 
% axis square; % grid on;
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',32,'Interpreter','Latex'); 
ylabel('concentration $C_i$ [kBq/mL]','FontSize',32,'Interpreter','Latex');
% title({'$C_i$'},'FontSize',25,'Interpreter','Latex')
fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'Ci_std'),'-dpng','-r300')