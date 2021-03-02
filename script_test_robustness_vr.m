clc
clear
close all

%% Supplementary code for testing the effect on BCM of different values of 
%% v=V_er/V_cyt ratio of the intracellular volumes of ER and cytosol

%% Define data path 
path_data = './data'; 
   % Data are available upon request to prof Gianmario Sambuceti (Sambuceti at unige.it).
path_results_suppl = './supplementary_results';
path_results = './results';
path_functions = './func';
addpath(path_functions)

mice_CT26 = {'#0 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#1 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#2 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#4 CT26 NO STS 2012.11.14 PRIMA PET';...
    '#1 CT26 NO STS 2012.11.20 SECONDA PET';...
    '#3 CT26 NO STS 2012.11.21 SECONDA PET'};

%% Define general parameters
% Number of repetition of the algorithm
tent = 50; 

% Value of v = V_er/V_cyt to be tested
v_all = 0.17*[1/2, 3/4, 1, 5/4, 2];
n_v = numel(v_all);

%% for each mouse
for i=2:length(mice_CT26)

disp(' '); disp(['mouse = ',mice_CT26{i}]); 
mouse = mice_CT26{i};

% Load data
if i <= 4
    folder = strcat('/PRIMA PET/',mouse,'/');
else
    folder = strcat('/SECONDA PET/',mouse,'/');
end  
[t,Ca,Ct] = import_voistat(path_data,folder);

% Load results from SKM to initialize k5
load(fullfile(path_results, sprintf('K_%s.mat', mouse)), 'K_Skf');
    
% Correct data
if i == 1 || i == 3
    Ca(21) = (Ca(20)+Ca(22))/2;
    Ct(21) = (Ct(20)+Ct(22))/2;
end
t = t'; Ca = @(tt)(interp1([0 t],[0 Ca'],tt,'linear',0));


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
    k5x = 1;
    while k5x > norm(K_Skf.mean)
        k5x = rand(1);
    end
    k6x = 0;
    
    [k1_BCM_vec(n, iv), k2_BCM_vec(n, iv), k3_BCM_vec(n, iv), k5_BCM_vec(n, iv), k6_BCM_vec(n, iv),...
        Cx_BCM_vec{n, iv}, Cxdata_BCM, relerr_BCM_vec(n, iv), iter_BCM_vec(n, iv)] = ...
        reconstruction_BCM(Ct,Ca,t,0,[0;0;0], k1x, k2x, k3x, k5x, k6x, v);

end

end

% store the results
K_BCM.k1 = k1_BCM_vec; K_BCM.k2 = k2_BCM_vec; K_BCM.k3 = k3_BCM_vec; K_BCM.k5 = k5_BCM_vec; K_BCM.k6 = k6_BCM_vec;
K_BCM.relerr = relerr_BCM_vec; K_BCM.iter = iter_BCM_vec; K_BCM.comp = Cx_BCM_vec;

%% Save
save(fullfile(path_results_suppl, sprintf('test_Vr_K_%s.mat', mouse)), ...
    'mouse','t','Ca','Ct','K_BCM');

end