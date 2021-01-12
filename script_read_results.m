clc
clear 
close all

path_results = '.';

mice_CT26 = {'#0 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#1 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#2 CT26 NO STS 2012.11.13 PRIMA PET';...
    '#4 CT26 NO STS 2012.11.14 PRIMA PET';...
    '#1 CT26 NO STS 2012.11.20 SECONDA PET';...
    '#3 CT26 NO STS 2012.11.21 SECONDA PET'};

n_mice = numel(mice_CT26);

k_name_bcm = {'k1', 'k2', 'k3', 'k5', 'k6'};
k_name_scm = {'k1', 'k2', 'k3', 'k4'};

for im = 1:n_mice
    
    mouse = mice_CT26{im};
    file_result = fullfile(path_results, strcat('K_',mouse,'.mat'));
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