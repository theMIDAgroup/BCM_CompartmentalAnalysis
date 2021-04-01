clc
clear
close all

%% Set path and general parameters
set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

path_figures = './figures';
path_functions = './func';
addpath(path_functions)

path_data = './data';

%% Define model parameters and load data
Vb = 0.15; % blood volume fraction (tumor)
Vi = 0.3; % interstizial volume fraction (tumor)
v = 0.17; % total volume of the endoplasmatic reticulum
Vr = v/(1+v);

alpha = [Vi+(1-Vr)*(1-Vb-Vi), (1-Vr)*(1-Vb-Vi), Vr*(1-Vb-Vi)];

t0 = 0;         % Inital data for the compartment model
C0 = [0; 0; 0];
y0 = [0; 0; 0];

load(fullfile(path_data, 'data_m1.mat'));
t = data.t'; Ca = data.Ca; clear data
Ca = @(tt)(interp1([0 t],[0 Ca'],tt,'linear',0));

%% Define reference values
k1_true = 0.32; k2_true = 0.37; k3_true = 0.45; k5_true = 0.51; k6_true = 0.03; 
    % Value estimated for mouse m1.
Mx = [[-(k2_true+k3_true);k3_true;0],[0;-k5_true;k5_true],[k6_true;0;-k6_true]];
Cx = concentration(k1_true,Mx,Ca,0,[0;0;0],t);
CT_true = (alpha*Cx + Vb*Ca(t))';

%% ********* Sensitivity analysis (one kinetic parameter at the time)
range = logspace(-1, 2, 20);

%% k1 
k1_vec = k1_true*range;
k2 = k2_true;
k3 = k3_true;
k5 = k5_true;
k6 = k6_true;

Cxdata = zeros(length(k1_vec), length(t));
for i = 1:length(k1_vec)
    k1 = k1_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(i,:) = (alpha*Cx + Vb*Ca(t))';
end

f_k1 = plot_3d_curves(t, k1_vec, Cxdata, [45, 15], 1, k1_true);
hold on
aux_t = [t'; t(end)];
plot3(aux_t, k1_true*ones(size(aux_t)), [CT_true; 0], 'k.--', 'linewidth', 4)
f_k1.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_3D_k1'),'-dpng','-r300')

%% k2
k1 = k1_true;
k2_vec = k2_true*range;
k3 = k3_true;
k5 = k5_true;
k6 = k6_true;

Cxdata = zeros(length(k2_vec), length(t));
for i = 1:length(k2_vec)
    k2 = k2_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(i,:) = (alpha*Cx + Vb*Ca(t))';
end

f_k2 = plot_3d_curves(t, k2_vec, Cxdata, [145, 15], 2, k2_true);
hold on
aux_t = [t'; t(end)];
plot3(aux_t, k2_true*ones(size(aux_t)), [CT_true; 0], 'k.--', 'linewidth', 4)
f_k2.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_3D_k2'),'-dpng','-r300')

%% k3
k1 = k1_true;
k2 = k2_true;
k3_vec = k3_true*range;
k5= k5_true;
k6 = k6_true;

Cxdata = zeros(length(k3_vec), length(t));
for i = 1:length(k3_vec)
    k3 = k3_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(i,:) = (alpha*Cx + Vb*Ca(t))';
end

f_k3 = plot_3d_curves(t, k3_vec, Cxdata, [45, 15], 3, k3_true);
hold on
aux_t = [t'; t(end)];
plot3(aux_t, k3_true*ones(size(aux_t)), [CT_true; 0], 'k.--', 'linewidth', 4)
f_k3.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_3D_k3'),'-dpng','-r300')

%% k5 
k1 = k1_true;
k2 = k2_true;
k3 = k3_true;
k5_vec = k5_true*range;
k6 = k6_true;

Cxdata = zeros(length(k5_vec), length(t));
for i = 1:length(k5_vec)
    k5 = k5_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(i,:) = (alpha*Cx + Vb*Ca(t))';
end

f_k5 = plot_3d_curves(t, k5_vec, Cxdata, [145, 15], 5, k5_true);
hold on
aux_t = [t'; t(end)];
plot3(aux_t, k5_true*ones(size(aux_t)), [CT_true; 0], 'k.--', 'linewidth', 4)
f_k5.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_3D_k5'),'-dpng','-r300')

%% k6
k1 = k1_true;
k2 = k2_true;
k3 = k3_true;
k5 = k5_true;
k6_vec = k6_true*range;

Cxdata = zeros(length(k6_vec), length(t));
for i = 1:length(k6_vec)
    k6 = k6_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(i,:) = (alpha*Cx + Vb*Ca(t))';
end

f_k6 = plot_3d_curves(t, k6_vec, Cxdata, [45, 15], 6, k6_true);
hold on
aux_t = [t'; t(end)];
plot3(aux_t, k6_true*ones(size(aux_t)), [CT_true; 0], 'k.--', 'linewidth', 4)
f_k6.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_3D_k6'),'-dpng','-r300')

%% ********* Relative local sensitivity S_T
% Initialization
range = logspace(-3, 3, 20);
n_scaling_factor = numel(range);
n_par = 5;

S_T = zeros(n_par, n_scaling_factor);

t = t';

%% k1
k2 = k2_true; k3 = k3_true; k5 = k5_true; k6 = k6_true;
for is = 1:n_scaling_factor
    k1 = k1_true*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t), [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k1 = f_comp_sol_activity(A_coeff, Ca(t), [1; 0; 0], t0, y0, t);
    ST_C_k1 = k1 * (y_k1 ./ C_time);
    
    CT_dk1 = alpha * y_k1;
    ST_k1 = abs(k1*(CT_dk1 ./ CT_time));
    
    S_T(1, is) = ST_k1(end);

end

clear k1 k2 k3 k5 k6

%% k2
k1 = k1_true; k3 = k3_true; k5 = k5_true; k6 = k6_true;
for is = 1:n_scaling_factor
    k2 = k2_true*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t), [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k2 = f_comp_sol_activity(A_coeff, -C_time(1, :)', [1; 0; 0], t0, y0, t);
    ST_C_k2 = abs(k2 * (y_k2 ./ C_time));
    
    CT_dk2 = alpha * y_k2;
    ST_k2 = abs(k2 * CT_dk2 ./ CT_time);
    
    S_T(2, is) = ST_k2(end);

end

clear k1 k2 k3 k5 k6

%% k3 
k1 = k1_true; k2 = k2_true; k5 = k5_true; k6 = k6_true;
for is = 1:n_scaling_factor
    k3 = k3_true*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t), [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k3 = f_comp_sol_activity(A_coeff, -C_time(1, :)', [1; 0; 0], t0, y0, t) + ...
        f_comp_sol_activity(A_coeff, C_time(1, :)', [0; 1; 0], t0, y0, t);
    ST_C_k3 = abs(k3 * (y_k3 ./ C_time));

    CT_dk3 = alpha * y_k3;
    ST_k3 = abs(k3 * CT_dk3 ./ CT_time);

    S_T(3, is) = ST_k3(end);

end

clear k1 k2 k3 k5 k6

%% k5
k1 = k1_true; k2 = k2_true; k3 = k3_true; k6 = k6_true;
for is = 1:n_scaling_factor
    k5 = k5_true*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t), [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k5 = f_comp_sol_activity(A_coeff, -C_time(2, :)', [0; 1; 0], t0, y0, t) + ...
        f_comp_sol_activity(A_coeff, C_time(2, :)', [0; 0; 1], t0, y0, t);
    ST_C_k5 = abs(k5 * (y_k5 ./ C_time));

    CT_dk5 = alpha * y_k5;
    ST_k5 = abs(k5 * CT_dk5 ./ CT_time);
    
    S_T(4, is) = ST_k5(end);

end

clear k1 k2 k3 k5 k6

%% k6
k1 = k1_true; k2 = k2_true; k3 = k3_true; k5 = k5_true;
for is = 1:n_scaling_factor
    k6 = k6_true*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t), [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k6 = f_comp_sol_activity(A_coeff, -C_time(3, :)', [1; 0; 0], t0, y0, t) + ...
        f_comp_sol_activity(A_coeff, C_time(3, :)', [0; 0; 1], t0, y0, t);
    ST_C_k6 = abs(k6 * (y_k6 ./ C_time));

    CT_dk6 = alpha * y_k6;
    ST_k6 = abs(k6 * CT_dk6 ./ CT_time);
    
    S_T(5, is) = ST_k6(end);

end

clear k1 k2 k3 k5 k6

%% Plot results
figure('units','normalized','outerposition',[0 0 0.5 0.7]);
set(gcf,'Color','none');
hold on
plot(range, S_T(1, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^1$')
plot(range, S_T(2, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^2$')
plot(range, S_T(3, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^3$')
plot(range, S_T(4, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^5$')
plot(range, S_T(5, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^6$')
lgd = legend('show', 'Location', 'Bestoutside');
lgd.FontSize = 25;
set(lgd, 'Interpreter', 'Latex');
set(gca, 'XScale', 'log', 'FontSize', 25, 'Xtick', [0.001, 0.01, 0.1, 1, 10, 100, 1000]);
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ylabel('$\overline{S}_T^j$', 'FontSize', 25, 'Interpreter', 'Latex')
xlabel('Scaling factor', 'FontSize', 25, 'Interpreter', 'Latex')
axis tight
print(fullfile(path_figures, 'ST_all_parameters'), '-dpng','-r300')

%% Auxiliary function for 3D plots
function f1 = plot_3d_curves(time, k_j_vec, Cxdata, view_par, idx_k, k_j_true)

    f1 = figure('units','normalized','outerposition',[0 0 0.5 0.7]);
    set(gcf,'Color','none');
    
    [aux_X, aux_Y] = meshgrid(time, k_j_vec);
    h = waterfall(aux_X, aux_Y, Cxdata);
    set(h, 'LineWidth', 4 );
    hidden off;
    view(view_par)
    zlim([0, 700])
    ylim([k_j_vec(1), k_j_vec(end)])
    aux_ytick = [k_j_vec(1), k_j_true, k_j_vec(end)];
    set(gca, 'yscale', 'log', 'ytick', aux_ytick)
    xlabel('Time [min]')
    ylabel(sprintf('$k_%d$', idx_k), 'Interpreter', 'latex')
    zlabel('concentration $C_T$ [kBq/mL]', 'Interpreter', 'latex')

    CD = get (h, 'CData');
    CD(end-2:end,:) = nan;
    set (h, 'CData', CD)
    
end




