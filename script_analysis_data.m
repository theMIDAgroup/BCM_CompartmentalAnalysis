clc
clear
close all

%% Main code to estimate the kinetic parameters of the BCM and the SCM.

%% Define data path 
path_data = './data'; 
path_functions = './func';
addpath(path_functions)
path_results = './results';

%% Initialization
% Number of mice
n_mice = 6;

% Number of repetition of the algorithm
tent = 50; 

% SCM
k1_Skf_vec = zeros(tent,1); k2_Skf_vec = zeros(tent,1);
k3_Skf_vec = zeros(tent,1); k4_Skf_vec = zeros(tent,1);
relerr_Skf_vec = zeros(tent,1); iter_Skf_vec = zeros(tent,1);
Cx_Skf_vec = cell(tent,1);

% BCM
k1_BCM_vec = zeros(tent,1); k2_BCM_vec = zeros(tent,1); 
k3_BCM_vec = zeros(tent,1); k5_BCM_vec = zeros(tent,1); k6_BCM_vec = zeros(tent,1);
relerr_BCM_vec = zeros(tent,1); iter_BCM_vec = zeros(tent,1);
Cx_BCM_vec = cell(tent,1);

%% for each mouse
for i=1:n_mice

fprintf(' Working with mouse m%d ', i)

% Load data
load(fullfile(path_data, sprintf('data_m%d', i)))
t = data.t; Ca = data.Ca; Ct = data.Ct; clear data

t = t'; Ca = @(tt)(interp1([0 t],[0 Ca'],tt,'linear',0));
    
%% SCM
disp(' '); disp('reconstruction Sokoloff...');
for n=1:tent
    
    disp(['n = ',num2str(n)]); 

    [k1_Skf_vec(n),k2_Skf_vec(n),k3_Skf_vec(n),k4_Skf_vec(n),...
        Cx_Skf_vec{n},Cxdata_Skf,relerr_Skf_vec(n),iter_Skf_vec(n)] = ...
        reconstruction_Skf(Ct,Ca,t,0,[0;0]);
    
end

% mean values of the reconstructed parameters
Km_Skf = [mean(k1_Skf_vec); mean(k2_Skf_vec); mean(k3_Skf_vec); mean(k4_Skf_vec)];
% standard deviations of the reconstructed parameters
Kstd_Skf = [std(k1_Skf_vec); std(k2_Skf_vec); std(k3_Skf_vec); std(k4_Skf_vec)];

% solve the direct problem to obtain the 'mean' compartment concentrations
M_Skf = [[-(Km_Skf(2)+Km_Skf(3));Km_Skf(3)],[Km_Skf(4);-Km_Skf(4)]];
Cxm_Skf = concentration(Km_Skf(1),M_Skf,Ca,0,[0;0],t);

% store the results
K_Skf.k1 = k1_Skf_vec; K_Skf.k2 = k2_Skf_vec; K_Skf.k3 = k3_Skf_vec; K_Skf.k4 = k4_Skf_vec;
K_Skf.relerr = relerr_Skf_vec; K_Skf.iter = iter_Skf_vec; K_Skf.comp = Cx_Skf_vec;
K_Skf.mean = Km_Skf; K_Skf.std = Kstd_Skf; K_Skf.comp_mean = Cxm_Skf;

%% BCM
disp(' '); disp('reconstruction BCM...');
for n=1:tent
    
    disp(['n = ',num2str(n)]); 
    
    % Initialize kinetic parameters
    k1x = rand(1);
    k2x = rand(1);
    k3x = rand(1);
    k5x = 1;
    while k5x > norm(K_Skf.mean)
        k5x = rand(1);
    end
    k6x = 0;
    
    [k1_BCM_vec(n),k2_BCM_vec(n),k3_BCM_vec(n),k5_BCM_vec(n),k6_BCM_vec(n),...
        Cx_BCM_vec{n},Cxdata_BCM,relerr_BCM_vec(n),iter_BCM_vec(n)] = ...
        reconstruction_BCM(Ct,Ca,t,0,[0;0;0], k1x, k2x, k3x, k5x, k6x);

end

% mean values of the reconstructed parameters
Km_BCM = [mean(k1_BCM_vec); mean(k2_BCM_vec); mean(k3_BCM_vec); mean(k5_BCM_vec); mean(k6_BCM_vec)];
% standard deviations of the reconstructed parameters
Kstd_BCM = [std(k1_BCM_vec); std(k2_BCM_vec); std(k3_BCM_vec); std(k5_BCM_vec); std(k6_BCM_vec)];

% solve the direct problem to obtain the 'mean' compartment concentrations
M_BCM = [[-(Km_BCM(2)+Km_BCM(3));Km_BCM(3);0],[0;-Km_BCM(4);Km_BCM(4)],[Km_BCM(5);0;-Km_BCM(5)]];
Cxm_BCM = concentration(Km_BCM(1),M_BCM,Ca,0,[0;0;0],t);

% store the results
K_BCM.k1 = k1_BCM_vec; K_BCM.k2 = k2_BCM_vec; K_BCM.k3 = k3_BCM_vec; K_BCM.k5 = k5_BCM_vec; K_BCM.k6 = k6_BCM_vec;
K_BCM.relerr = relerr_BCM_vec; K_BCM.iter = iter_BCM_vec; K_BCM.comp = Cx_BCM_vec;
K_BCM.mean = Km_BCM; K_BCM.std = Kstd_BCM; K_BCM.comp_mean = Cxm_BCM;

%% Save
save(fullfile(path_results, sprintf('K_m%d.mat', i)), ...
            't','Ca','Ct','K_Skf','K_BCM');

end












