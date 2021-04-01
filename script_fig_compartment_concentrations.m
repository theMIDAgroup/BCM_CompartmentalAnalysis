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

load(fullfile(path_results, 'K_m1.mat'));

% Colors
dark_blue = [0 0 0.8];
orange = [1 0.5 0];
gold = [1 0.8 0];
dark_green = [0 0.8 0];

light_red = [0.8 0 0];
light_blue =[0, 250, 256] / 256;
light_green = [157,204,0] / 256;
light_orange = [153,0,0]/256;

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
CT_BCM_std = (1/sqrt(n_run))*std(temp_CT_BCM, [], 2);


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
CT_SKM_std = (1/sqrt(n_run))*std(temp_CT_SKM, [], 2);

%% Compute error in reconstructing C_T
rel_err_BCM = vecnorm(temp_CT_BCM-repmat(Ct, 1, n_run), 2, 1) / norm(Ct);
rel_err_SKM = vecnorm(temp_CT_SKM-repmat(Ct, 1, n_run), 2, 1) / norm(Ct);

fprintf('Relative error BCM = %2.3f +/- %2.3f (max = %2.3f) \n', ...
    mean(rel_err_BCM), std(rel_err_BCM), max(rel_err_BCM))
fprintf('Relative error SCM = %2.3f +/- %2.3f (max = %2.3f) \n', ...
    mean(rel_err_SKM), std(rel_err_SKM), max(rel_err_SKM))

%% Plots
%% P.a. BCM fit of the total concentration CT
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 4, 'Marker','o', 'Markersize',5);
hold on
errorbar(t, CT_BCM_mean, CT_BCM_std, 'Color', 'k', 'LineStyle','-', ...
    'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
axis([0 40 0 350]); 
axis square;
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Experimental data', 'Reconstructed'},'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'BCM_CT'),'-dpng','-r300')

%% P.b. SKM fit of the total concentration CT
figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
plot(t, Ct,'Color', dark_green, 'Linewidth', 4, 'Marker','o', 'Markersize',5);
hold on
errorbar(t, CT_SKM_mean, CT_SKM_std, 'Color', 'k', 'LineStyle','-', ...
    'Linewidth', 2.5, 'Marker','o', 'Markersize',5);
axis([0 40 0 350]); 
axis square;
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Experimental data', 'Reconstructed'},'FontSize',25,'Location','southeast');
print(fullfile(path_figures, 'SCM_CT'),'-dpng','-r300')

%% Pc. BCM compartments
C_BCM_mean = mean(temp_C_BCM, 3);
C_BCM_std = (1/sqrt(n_run))*std(temp_C_BCM, [], 3);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
errorbar(t, C_BCM_mean(1, :), C_BCM_std(1, :), 'Color',dark_blue, 'LineWidth',3);
hold on;
errorbar(t,  C_BCM_mean(2, :), C_BCM_std(2, :), 'Color',orange, 'LineWidth',3);
errorbar(t,  C_BCM_mean(3, :), C_BCM_std(3, :), 'Color',gold,'LineWidth',3);
axis([0 40 0 1500]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$C_f$',...
    '$C_p$','$C_r$'}, ...
    'FontSize',25,'Location','northwest');
print(fullfile(path_figures, 'BCM_comp'),'-dpng','-r300')

%% P2.e. SCM compartments
C_SKM_mean = mean(temp_C_SKM, 3);
C_SKM_std = (1/sqrt(n_run))*std(temp_C_SKM, [], 3);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
errorbar(t,  C_SKM_mean(1, :), C_SKM_std(1, :), 'Color',dark_blue,'LineWidth',3);
hold on;
errorbar(t, C_SKM_mean(2, :),  C_SKM_std(2, :), 'Color',orange,'LineWidth', 3);
axis([0 40 0 400]); 
axis square;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'$C_f$','$C_p$'},...
    'FontSize',25,'Location','northwest');
fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'SCM_comp'),'-dpng','-r300')

%% P2.f. BCM activity in different volumes
act_blood = Vb*Ca(t);
act_int = squeeze(Vi*temp_C_BCM(1, :, :));
act_cyt = squeeze(alpha_BCM(2)*temp_C_BCM(1, :, :) + alpha_BCM(2)*temp_C_BCM(2, :, :));
act_er = squeeze(alpha_BCM(3)*temp_C_BCM(3, :, :));

Ca_minus = act_blood - 2*sqrt(act_blood); 
act_blood_std = abs(act_blood - Ca_minus);
act_int_mean = mean(act_int, 2);
act_int_std = (1/sqrt(n_run))*std(act_int, [], 2);
act_cyt_mean = mean(act_cyt, 2);
act_cyt_std = (1/sqrt(n_run))*std(act_cyt, [], 2);
act_er_mean = mean(act_er, 2);
act_er_std = (1/sqrt(n_run))*std(act_er, [], 2);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
errorbar(t, act_int_mean, act_int_std, 'Color',light_blue,...
    'LineWidth',3,'Markersize',13);
hold on
errorbar(t,  act_cyt_mean, act_cyt_std, 'Color',light_orange,...
    'LineWidth',3,'Markersize',13);
errorbar(t,  act_er_mean, act_er_std, 'Color', gold,...
    'LineWidth',3,'Markersize',13);
axis([0 40 0 250]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('weighted concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Interstitium', 'Cytosol', 'Reticulum '}, ...
    'FontSize',25,'Location','northwest');
print(fullfile(path_figures, 'BCM_volumes'),'-dpng','-r300')

%% P2.f. SKM activity in different volumes
act_int_skm = squeeze(Vi*temp_C_SKM(1, :, :));
act_cyt_skm = squeeze(alpha_SKM(2)*temp_C_SKM(1, :, :)+alpha_SKM(2)*temp_C_SKM(2, :, :));

act_int_skm_mean = mean(act_int_skm, 2);
act_int_skm_std = (1/sqrt(n_run))*std(act_int_skm, [], 2);

act_cyt_skm_mean = mean(act_cyt_skm, 2);
act_cyt_skm_std = (1/sqrt(n_run))*std(act_cyt_skm, [], 2);

figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf, 'Color','none');
errorbar(t, act_int_skm_mean, act_int_skm_std, 'Color',light_blue,...
    'LineWidth',3,'Markersize',13);
hold on
errorbar(t,  act_cyt_skm_mean, act_cyt_skm_std, 'Color',light_orange,...
    'LineWidth',3,'Markersize',13);
axis([0 40 0 250]); 
axis square; %grid on;
set(gca,'xtick',0:10:40);
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('weighted concentration [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend({'Interstistium', 'Cytosol'}, ...
    'FontSize',25,'Location','northwest');
print(fullfile(path_figures, 'SKM_volumes'),'-dpng','-r300')









           




