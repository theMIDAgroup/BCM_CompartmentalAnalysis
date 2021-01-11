clc
clear
close all

set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

%% Define paths and load data
path_results = './results/';
path_figures = './figures/';

load(fullfile(path_results, 'K_#0 CT26 NO STS 2012.11.13 PRIMA PET.mat'));

% Colors
dark_blue = [0 0 0.8];
orange = [1 0.5 0];
gold = [1 0.8 0];

%% BCM compartments 
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t,K_BCM.comp_mean(1,:),'Color',dark_blue,...
    'Marker','x','LineWidth',3,'Markersize',13);
hold on;
plot(t,K_BCM.comp_mean(2,:),'Color',orange,...
    'Marker','*','LineWidth',3,'Markersize',13);
plot(t,K_BCM.comp_mean(3,:),'Color',gold,...
    'Marker','s','LineWidth',3,'Markersize',13);
axis([0 40 0 1400]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$C_f$','$C_p$','$C_r$'},'FontSize',25,'Location','northeastoutside');
print(fullfile(path_figures, 'BCM_comp_tissue'),'-dpng','-r300')

%% SCM compartments
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t,K_Skf.comp_mean(1,:),'Color',dark_blue,...
    'Marker','x','LineWidth',3,'Markersize',13);
hold on;
plot(t,K_Skf.comp_mean(2,:),'Color',orange,...
    'Marker','*','LineWidth',3,'Markersize',13);
axis square;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$C_f$','$C_p$'},'FontSize',25,'Location','northeastoutside');
fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'SCM_comp_tissue'),'-dpng','-r300')
