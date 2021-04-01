clc
clear
close all

set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

%% Define data path 
path_results = './results/';
path_functions = './func';
addpath(path_functions)

%% Define general parameters
Vb = 0.15; Vi = 0.3;
t0 = 0;
C0_BCM = [0; 0; 0];

im = 1; % Chose mouse

%   Load data
load(fullfile(path_results, sprintf('test_Vr_K_m%d.mat', im)))
n_run = size(K_BCM.k1, 1);
n_time_point = size(t, 2);

v_all = K_BCM.v_all;
vr_all = v_all ./ (1+v_all);

%% Compute mean and standard deviation of the kinetic parameters
k1_mean = zeros(numel(v_all), 1); k1_std = k1_mean;
k2_mean = k1_mean; k2_std = k1_mean;
k3_mean = k1_mean; k3_std = k1_mean;
k5_mean = k1_mean; k5_std = k1_mean;
k6_mean = k1_mean; k6_std = k1_mean;

for iv = 1:numel(v_all)
   k1_mean(iv) = mean(K_BCM.k1(:, iv), 1);
   k1_std(iv) = std(K_BCM.k1(:, iv), [], 1);

   k2_mean(iv) = mean(K_BCM.k2(:, iv), 1);
   k2_std(iv) = std(K_BCM.k2(:, iv), [], 1);

   k3_mean(iv) = mean(K_BCM.k3(:, iv), 1);
   k3_std(iv) = std(K_BCM.k3(:, iv), [], 1);

   k5_mean(iv) = mean(K_BCM.k5(:, iv), 1);
   k5_std(iv) = std(K_BCM.k5(:, iv), [], 1);

   k6_mean(iv) = mean(K_BCM.k6(:, iv), 1);
   k6_std(iv) = std(K_BCM.k6(:, iv), [], 1);
end

%% Table for the kinetick parameters
path_table = fullfile(path_results, sprintf('table_k_m%d.txt', im));
fileID = fopen(path_table, 'w');
for iv = 1:numel(v_all)
    fprintf(fileID, '%1.2f & %1.2f & %1.2f $\\pm$ %1.3f  & %1.2f $\\pm$ %1.2f & %1.2f $\\pm$ %1.2f & %1.2f $\\pm$ %1.2f & %1.2f $\\pm$ %1.2f \\\\ \\hline \n', ...
        v_all(iv), vr_all(iv), k1_mean(iv), k1_std(iv), k2_mean(iv), k2_std(iv), k3_mean(iv), k3_std(iv), ...
        k5_mean(iv), k5_std(iv), k6_mean(iv), k6_std(iv));
end
fclose(fileID);

