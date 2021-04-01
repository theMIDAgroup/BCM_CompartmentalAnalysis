clc
clear
close all

%% Supplementary code for testing the effect on BCM of an error in the value 
%% of v=V_er/V_cyt (ratio of the intracellular volumes of ER and cytosol)

% Note: this code needs the results from the SCM provided by the code 
%       'script_analysis_data.m'.

%% Define data path 
path_data = './data'; 
path_results = './results';

path_functions = './func';
addpath(path_functions)

%% Define general parameters
% Number of mice
n_mice = 6;

% Number of repetition of the algorithm
tent = 50; 

% Value of v = V_er/V_cyt to be tested
v_all = 0.17*[1/2, 3/4, 1, 5/4, 3/2];
n_v = numel(v_all);

%% for each mouse
for i=1:n_mice

fprintf(' Working with mouse m%d ', i)

% Load data
load(fullfile(path_data, sprintf('data_m%d', i)))
t = data.t; Ca = data.Ca; Ct = data.Ct; clear data

t = t'; Ca = @(tt)(interp1([0 t],[0 Ca'],tt,'linear',0));

% Load results from the SCM to initialize k5
load(fullfile(path_results, sprintf('K_m%d', i)), 'K_Skf')

%% BCM
% Intialize output variables
k1_BCM_vec = zeros(tent, n_v); k2_BCM_vec = zeros(tent, n_v); 
k3_BCM_vec = zeros(tent, n_v); k5_BCM_vec = zeros(tent, n_v);
k6_BCM_vec = zeros(tent, n_v);
relerr_BCM_vec = zeros(tent, n_v); iter_BCM_vec = zeros(tent, n_v);
Cx_BCM_vec = cell(tent, n_v);

for iv = 1:n_v

v = v_all(iv);
disp(' '); fprintf('Testing v = %2.2f \n', v);

for n=1:tent
    
    disp(['n = ',num2str(n)]); 
    
    % Initialize kinetic parameters
    k1x = rand(1);
    k2x = rand(1);
    k3x = rand(1);
    k5x = rand(1)* norm(K_Skf.mean);
    k6x = 0;
    
    [k1_BCM_vec(n, iv), k2_BCM_vec(n, iv), k3_BCM_vec(n, iv), k5_BCM_vec(n, iv), k6_BCM_vec(n, iv),...
        Cx_BCM_vec{n, iv}, Cxdata_BCM, relerr_BCM_vec(n, iv), iter_BCM_vec(n, iv)] = ...
        reconstruction_BCM(Ct,Ca,t,0,[0;0;0], k1x, k2x, k3x, k5x, k6x, v);

end

end

% store the results
K_BCM.k1 = k1_BCM_vec; K_BCM.k2 = k2_BCM_vec; K_BCM.k3 = k3_BCM_vec; K_BCM.k5 = k5_BCM_vec; K_BCM.k6 = k6_BCM_vec;
K_BCM.relerr = relerr_BCM_vec; K_BCM.iter = iter_BCM_vec; K_BCM.comp = Cx_BCM_vec; ...
K_BCM.v_all = v_all;

%% Save
save(fullfile(path_results, sprintf('test_Vr_K_m%d.mat', i)), ...
             't','Ca','Ct','K_BCM');

end