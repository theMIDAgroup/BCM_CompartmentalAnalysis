clc
clear
close all

set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');


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
dark_green = [0 0.8 0];

%% Define general parameters
t0 = 0;
C0_BCM = [0; 0; 0];
C0_SKM = [0; 0];

Vb = 0.15; Vi = 0.3;
v = 0.17; Vr = v/(1+v);
alpha_BCM = [Vi+(1-Vr)*(1-Vb-Vi),(1-Vr)*(1-Vb-Vi),Vr*(1-Vb-Vi)];
alpha_SKM = [1-Vb,1-Vb-Vi];

n_run = size(K_BCM.k1, 1);
n_time_point = size(t, 2);

%% Compute reconstructed CT (from mean values of k)
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


%% Compute reconstructed CT (average of the curves at each run)
% - BCM (all run)
temp_C_BCM = zeros(3, n_time_point, n_run);
temp_CT_BCM = zeros(n_time_point, n_run);
for ir = 1:n_run
    k1_BCM = K_BCM.k1(ir); 
    k2_BCM = K_BCM.k2(ir);
    k3_BCM = K_BCM.k3(ir);
    k5_BCM = K_BCM.k5(ir);
    k6_BCM = K_BCM.k6(ir);
    
    M_BCM = [[-(k2_BCM+k3_BCM);k3_BCM;0],[0;-k5_BCM;k5_BCM],[k6_BCM;0;-k6_BCM]];
    temp_C_BCM(:, :, ir) = concentration(k1_BCM, M_BCM, Ca, t0, C0_BCM, t);
    temp_CT_BCM(:, ir) = (alpha_BCM*temp_C_BCM(:, :, ir)+ Vb*Ca(t))';
end

CT_BCM_mean = mean(temp_CT_BCM, 2);
CT_BCM_std = std(temp_CT_BCM, [], 2);


% - SKM (all run)
temp_C_SKM = zeros(2, n_time_point, n_run);
temp_CT_SKM = zeros(n_time_point, n_run);
for ir = 1:n_run
    k1_SKM = K_Skf.k1(ir);
    k2_SKM = K_Skf.k2(ir);
    k3_SKM = K_Skf.k3(ir);
    k4_SKM = K_Skf.k4(ir);
    
    M_SKM = [[-(k2_SKM+k3_SKM);k3_SKM],[k4_SKM;-k4_SKM]];
    temp_C_SKM(:, :, ir) = concentration(k1_SKM, M_SKM, Ca, t0, C0_SKM, t);
    temp_CT_SKM(:, ir) = (alpha_SKM*temp_C_SKM(:, :, ir) + Vb*Ca(t))';
end
CT_SKM_mean = mean(temp_CT_SKM, 2);
CT_SKM_std = std(temp_CT_SKM, [], 2);

%% ***************** Plot 1. Mean values of k *****************************

%% P1.a. BCM fit of the total concentration CT256
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
plot(t, CT_BCM, 'Color', 'k', 'Linewidth', 2.5, 'Marker','o', 'Markersize',5);

axis square; %grid on;
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'BCM_CT_avek'),'-dpng','-r300')

%% P1.b. SKF fit of the total concentration CT
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
plot(t, CT_SKM, 'Color', 'k', 'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
axis square; %grid on;
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'SCM_CT_avek'),'-dpng','-r300')

%% P1.c. BCM compartments
temp_C = alpha_BCM' .* C_BCM;
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, temp_C(1,:),'Color',dark_blue,...
    'Marker','x','LineWidth',3,'Markersize',13);
hold on;
plot(t, temp_C(2,:),'Color',orange,...
    'Marker','*','LineWidth',3,'Markersize',13);
plot(t, temp_C(3,:),'Color',gold,...
    'Marker','s','LineWidth',3,'Markersize',13);
axis([0 40 0 200]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$ \frac{V_{int} + V_{cyt}}{V_{tot}} C_f$',...
    '$ \frac{V_{cyt}}{V_{tot}}C_p$','$\frac{V_{er}}{V_{tot}}C_r$'}, ...
    'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'BCM_comp_tissue_avek'),'-dpng','-r300')

%% P1.d. SKM compartments
temp_C = alpha_SKM' .* C_SKM;
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, temp_C(1,:),'Color',dark_blue,...
    'Marker','x','LineWidth',3,'Markersize',13);
hold on;
plot(t, temp_C(2,:),'Color',orange,...
    'Marker','*','LineWidth',3,'Markersize',13);
axis([0 40 0 200]); 
axis square;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$\frac{V_{int} + V_{cyt}}{V_{tot}} C_f$','$ \frac{V_{cyt}}{V_{tot}} C_p$'},...
    'FontSize',25,'Location','southeast');
fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'SCM_comp_tissue_avek'),'-dpng','-r300')

%% **************** Plot 2. Average of the curves *************************
%% P2.a. BCM fit of the total concentration CT
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
errorbar(t, CT_BCM_mean, CT_BCM_std, 'Color', 'k', 'LineStyle','-', ...
    'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
hold on
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'BCM_CT_avecurves'),'-dpng','-r300')

%% P2.b. SKM fit of the total concentration CT
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 3.5, 'Marker','o', 'Markersize',5);
hold on
errorbar(t, CT_SKM_mean, CT_SKM_std, 'Color', 'k', 'LineStyle','-', ...
    'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Measured', 'Reconstructed'},'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'SCM_CT_avecurves'),'-dpng','-r300')

%% P2.c. BCM compartments
C_BCM_mean = mean(alpha_BCM' .* temp_C_BCM, 3);
C_BCM_std = std(alpha_BCM' .* temp_C_BCM, [], 3);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
errorbar(t, C_BCM_mean(1, :), C_BCM_std(1, :), 'Color',dark_blue,...
    'Marker','x','LineWidth',3,'Markersize',13);
hold on;
errorbar(t,  C_BCM_mean(2, :), C_BCM_std(2, :), 'Color',orange,...
    'Marker','*','LineWidth',3,'Markersize',13);
errorbar(t,  C_BCM_mean(3, :), C_BCM_std(3, :), 'Color',gold,...
    'Marker','s','LineWidth',3,'Markersize',13);
axis([0 40 0 200]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$ \frac{V_{int} + V_{cyt}}{V_{tot}} C_f$',...
    '$ \frac{V_{cyt}}{V_{tot}}C_p$','$\frac{V_{er}}{V_{tot}}C_r$'}, ...
    'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'BCM_comp_tissue_avecurves'),'-dpng','-r300')

%% P2.e. SCM compartments
C_SKM_mean = mean(alpha_SKM' .* temp_C_SKM, 3);
C_SKM_std = std(alpha_SKM' .* temp_C_SKM, [], 3);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
errorbar(t,  C_SKM_mean(1, :), C_SKM_std(1, :), 'Color',dark_blue,...
    'Marker','x','LineWidth',3,'Markersize',13);
hold on;
errorbar(t, C_SKM_mean(2, :),  C_SKM_std(2, :), 'Color',orange,...
    'Marker','*','LineWidth',3,'Markersize',13);
axis([0 40 0 200]); 
axis square;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$\frac{V_{int} + V_{cyt}}{V_{tot}} C_f$','$ \frac{V_{cyt}}{V_{tot}} C_p$'},...
    'FontSize',25,'Location','southeast');
fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'SCM_comp_tissue_avecurves'),'-dpng','-r300')

%% P2.f. BCM activity in different volumes
act_blood = Vb*Ca(t);
act_int = squeeze(Vi*temp_C_BCM(1, :, :));
act_cyt = squeeze(alpha_BCM(2)*temp_C_BCM(1, :, :) + alpha_BCM(2)*temp_C_BCM(2, :, :));
act_er = squeeze(alpha_BCM(3)*temp_C_BCM(3, :, :));

act_int_mean = mean(act_int, 2);
act_int_std = std(act_int, [], 2);
act_cyt_mean = mean(act_cyt, 2);
act_cyt_std = std(act_cyt, [], 2);
act_er_mean = mean(act_er, 2);
act_er_std = std(act_er, [], 2);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, act_blood, 'Color', dark_green, 'Linewidth', 2)
hold on
errorbar(t, act_int_mean, act_int_std, 'Color',dark_blue,...
    'LineWidth',3,'Markersize',13);
errorbar(t,  act_cyt_mean, act_cyt_std, 'Color',orange,...
    'LineWidth',3,'Markersize',13);
errorbar(t,  act_er_mean, act_er_std, 'Color',gold,...
    'Marker','s','LineWidth',3,'Markersize',13);
axis([0 40 0 260]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'blood', 'Interstistial', 'cytosol', 'reticulum'}, ...
    'FontSize',25,'Location','northwest');
print(fullfile(path_figures, 'BCM_volumes_avecurves'),'-dpng','-r300')

%% P2.f. SKM activity in different volumes
act_int_skm = squeeze(Vi*temp_C_SKM(1, :, :));
act_cyt_skm = squeeze(alpha_SKM(2)*temp_C_SKM(1, :, :)+alpha_SKM(2)*temp_C_SKM(2, :, :));

act_int_skm_mean = mean(act_int_skm, 2);
act_int_skm_std = std(act_int_skm, [], 2);

act_cyt_skm_mean = mean(act_cyt_skm, 2);
act_cyt_skm_std = std(act_cyt_skm, [], 2);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, act_blood, 'Color', dark_green, 'Linewidth', 2)
hold on
errorbar(t, act_int_skm_mean, act_int_skm_std, 'Color',dark_blue,...
    'LineWidth',3,'Markersize',13);
errorbar(t,  act_cyt_skm_mean, act_cyt_skm_std, 'Color',orange,...
    'LineWidth',3,'Markersize',13);
axis([0 40 0 260]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'blood', 'Interstistial', 'cytosol', 'reticulum'}, ...
    'FontSize',25,'Location','northwest');
print(fullfile(path_figures, 'SKM_volumes_avecurves'),'-dpng','-r300')









           




