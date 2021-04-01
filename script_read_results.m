clc
clear 
close all

%% Read results and generate Table 2 and Table 3 (estimated values of the
%  kinetic parameters.

path_results = './results';

n_mice =6;

k_name_bcm = {'k1', 'k2', 'k3', 'k5', 'k6'};
k_name_scm = {'k1', 'k2', 'k3', 'k4'};

for im = 1:n_mice
    
    file_result = fullfile(path_results, sprintf('K_m%d.mat', im));
    load(file_result)
    
    fprintf('******* MOUSE m%d *******  \n', im)
    fprintf('BCM \n')
    for ik = 1:numel(k_name_bcm)
        fprintf('%s - mean: %2.2f std = %2.2f \n', ...
            k_name_bcm{ik}, K_BCM.mean(ik), K_BCM.std(ik))
    end
    
    fprintf('SCM \n')
    for ik = 1:numel(k_name_scm)
        fprintf('%s - mean: %2.2f std = %2.2f \n', ...
            k_name_scm{ik}, K_Skf.mean(ik), K_Skf.std(ik))
    end
    
    
end