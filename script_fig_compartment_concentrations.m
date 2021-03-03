clc
clear
close all

set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

%% Define color
dark_green = [0 0.8 0];

%% Define paths and load data
path_results = './results/';
path_figures = './figures/';
path_functions = './func';
addpath(path_functions)

load(fullfile(path_results, 'K_#0 CT26 NO STS 2012.11.13 PRIMA PET.mat'));

% Colors
dark_blue = [0 0 0.8];
orange = [1 0.5 0];
gold = [1 0.8 0];

%% Compute reconstructed CT

% - General parameters
t0 = 0;
C0_BCM = [0; 0; 0];
C0_SKM = [0; 0];

Vb = 0.15; Vi = 0.3;
v = 0.17; Vr = v/(1+v);
alpha_BCM = [Vi+(1-Vr)*(1-Vb-Vi),(1-Vr)*(1-Vb-Vi),Vr*(1-Vb-Vi)];
alpha_SKM = [1-Vb,1-Vb-Vi];

n_run = size(K_BCM.k1, 1);
n_time_point = size(t, 2);

% - BCM (table 2)
k1_BCM = K_BCM.mean(1); 
k2_BCM = K_BCM.mean(2);
k3_BCM = K_BCM.mean(3);
k5_BCM = K_BCM.mean(4);
k6_BCM = K_BCM.mean(5);

M_BCM = [[-(k2_BCM+k3_BCM);k3_BCM;0],[0;-k5_BCM;k5_BCM],[k6_BCM;0;-k6_BCM]];
C_BCM = concentration(k1_BCM, M_BCM, Ca, t0, C0_BCM, t);
CT_BCM = (alpha_BCM*C_BCM + Vb*Ca(t))';

% - SKM (table 3)
k1_SKM = K_Skf.mean(1);
k2_SKM = K_Skf.mean(2);
k3_SKM = K_Skf.mean(3);
k4_SKM = K_Skf.mean(4);

M_SKM = [[-(k2_SKM+k3_SKM);k3_SKM],[k4_SKM;-k4_SKM]];
C_SKM = concentration(k1_SKM, M_SKM, Ca, t0, C0_SKM, t);
CT_SKM = (alpha_SKM*C_SKM + Vb*Ca(t))';

% - BCM all run 
temp_CT_BCM = zeros(n_time_point, n_run);
for ir = 1:n_run
    k1_BCM = K_BCM.k1(ir); 
    k2_BCM = K_BCM.k2(ir);
    k3_BCM = K_BCM.k3(ir);
    k5_BCM = K_BCM.k5(ir);
    k6_BCM = K_BCM.k6(ir);
    
    M_BCM = [[-(k2_BCM+k3_BCM);k3_BCM;0],[0;-k5_BCM;k5_BCM],[k6_BCM;0;-k6_BCM]];
    C_BCM = concentration(k1_BCM, M_BCM, Ca, t0, C0_BCM, t);
    temp_CT_BCM(:, ir) = (alpha_BCM*C_BCM + Vb*Ca(t))';
end
CT_BCM_mean = mean(temp_CT_BCM, 2);
CT_BCM_std = std(temp_CT_BCM, [], 2);

% - SKM all run
temp_CT_SKM = zeros(n_time_point, n_run);
for ir = 1:n_run
    k1_SKM = K_Skf.k1(ir);
    k2_SKM = K_Skf.k2(ir);
    k3_SKM = K_Skf.k3(ir);
    k4_SKM = K_Skf.k4(ir);
    
    M_SKM = [[-(k2_SKM+k3_SKM);k3_SKM],[k4_SKM;-k4_SKM]];
    C_SKM = concentration(k1_SKM, M_SKM, Ca, t0, C0_SKM, t);
    temp_CT_SKM(:, ir) = (alpha_SKM*C_SKM + Vb*Ca(t))';
end
CT_SKM_mean = mean(temp_CT_SKM, 2);
CT_SKM_std = std(temp_CT_SKM, [], 2);

%% BCM fit of the total concentration CT
%   - Mean values of the kinetic parameters
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
% set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
plot(t, CT_BCM, 'Color', 'k', 'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
axis square; %grid on;
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');

%   - Mean value of the curves
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
% set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
errorbar(t, CT_BCM_mean, CT_BCM_std, 'Color', 'k', 'LineStyle','-', ...
    'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
hold on
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');


%% SKM fit of the total concentration CT
%   - Mean values of the kinetic parameters
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
% set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
plot(t, CT_SKM, 'Color', 'k', 'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
axis square; %grid on;
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');

%   - Mean values of the curves
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
% set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
errorbar(t, CT_SKM_mean, CT_SKM_std, 'Color', 'k', 'LineStyle','-', ...
    'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');

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
