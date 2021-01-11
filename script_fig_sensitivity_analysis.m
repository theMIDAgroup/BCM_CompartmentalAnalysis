clc
clear
close all

set(0,'DefaultAxesFontSize',25)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

path_figures = './figures';
path_functions = './func';
addpath(path_functions)

%% Define model parameters
Vb = 0.15; % blood volume fraction (tumor)
Vi = 0.3; % interstizial volume fraction (tumor)
v = 0.17; % total volume of the enoplasmatic reticulum
Vr = v/(1+v); 

alpha = [Vi+(1-Vr)*(1-Vb-Vi), (1-Vr)*(1-Vb-Vi), Vr*(1-Vb-Vi)];

% Acquisition time
time = [7.5 22.50262 37.50008 52.50036 67.50043 82.50053 97.50076 ...
    112.5008 127.501005 142.501 165.0011 195.0016 225.0018 255.002 ...
    285.0021 334.1693 394.1695 454.1707 514.171015 574.1711 685.1639 ...
    835.1641 1065.2923 1370.788 1675.5025 1980.81945 2285.738];
t = (time./60)';
ntime = length(t);

% simulated Input Function
Ca = @(tt)(6.6991e+06*(tt.^5.4135e+00).*exp(-tt/9.65e-02)+3.611738e+02*exp(-tt/2.85452e+01)).';

%% Sensitivity analysis (one kinetic parameter at the time)
range = [0.1, 0.5, 1, 2, 10, 20, 100];

%% k1
k1_vec = 0.3*range;
k2 = 0.4;
k3 = 0.6;
k5= 0.6;
k6 = 0.06;

Cxdata = zeros(length(t),length(k1_vec));
for i = 1:length(k1_vec)
    k1 = k1_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(:,i) = (alpha*Cx + Vb*Ca(t))';
end

figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf,'Color','none');
hold on;
axis square;% grid on;
for i = 1:length(k1_vec)    
    if range(i) == 1
        plot(t,Cxdata(:,i), '-x', 'Markersize', 10, 'LineWidth', 3.5, ...
            'Displayname', strcat(['k1 = ', num2str(k1_vec(i))]));
    else
    plot(t,Cxdata(:,i), 'LineWidth', 3.5, 'Displayname', strcat(['k1 = ', num2str(k1_vec(i))]));   
    end
end
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend('show', 'Location','northeastoutside')
ylim([0, 500])
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_k1'),'-dpng','-r300')

%% k2
k1 = 0.3;
k2_vec = 0.4*range;
k3 = 0.6;
k5= 0.6;
k6 = 0.06;

Cxdata = zeros(length(t),length(k2_vec));
for i = 1:length(k2_vec)
    k2 = k2_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(:,i) = (alpha*Cx + Vb*Ca(t))';
end

figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf,'Color','none');
hold on;
axis square;% grid on;
for i = 1:length(k2_vec)    
    if range(i) == 1
        plot(t,Cxdata(:,i), '-x', 'Markersize', 10, 'LineWidth', 3.5, ...
            'Displayname', strcat(['k2 = ', num2str(k2_vec(i))]));
    else
    plot(t,Cxdata(:,i), 'LineWidth', 3.5, 'Displayname', strcat(['k2 = ', num2str(k2_vec(i))]));   
    end
end
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend('show', 'Location','northeastoutside')
ylim([0, 500])
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_k2'),'-dpng','-r300')

%% k3
k1 = 0.3;
k2 = 0.4;
k3_vec = 0.6*range;
k5= 0.6;
k6 = 0.06;

Cxdata = zeros(length(t),length(k3_vec));
for i = 1:length(k3_vec)
    k3 = k3_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(:,i) = (alpha*Cx + Vb*Ca(t))';
end

figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf,'Color','none');
hold on;
axis square;% grid on;
for i = 1:length(k3_vec)    
    if range(i) == 1
        plot(t,Cxdata(:,i), '-x', 'Markersize', 10, 'LineWidth', 3.5, ...
            'Displayname', strcat(['k3 = ', num2str(k3_vec(i))]));
    else
    plot(t,Cxdata(:,i), 'LineWidth', 3.5, 'Displayname', strcat(['k3 = ', num2str(k3_vec(i))]));   
    end
end
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend('show', 'Location','northeastoutside')
ylim([0, 500])
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_k3'),'-dpng','-r300')

%% k5 
k1 = 0.3;
k2 = 0.4;
k3 = 0.6;
k5_vec = 0.6*range;
k6 = 0.06;

Cxdata = zeros(length(t),length(k5_vec));
for i = 1:length(k5_vec)
    k5 = k5_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(:,i) = (alpha*Cx + Vb*Ca(t))';
end

figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf,'Color','none');
hold on;
axis square;% grid on;
for i = 1:length(k5_vec)    
       
    if range(i) == 1
        plot(t,Cxdata(:,i), '-x', 'Markersize', 10, 'LineWidth', 3.5, ...
            'Displayname', strcat(['k5 = ', num2str(k5_vec(i))]));
    else
        plot(t,Cxdata(:,i), 'LineWidth', 3.5, 'Displayname', strcat(['k5 = ', num2str(k5_vec(i))]));
    end
end
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend('show', 'Location','northeastoutside')
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ylim([0, 500])

fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_k5'),'-dpng','-r300')

%% k6
k1 = 0.3;
k2 = 0.4;
k3 = 0.6;
k5= 0.6;
k6_vec = 0.06*range;

Cxdata = zeros(length(t),length(k6_vec));
for i = 1:length(k6_vec)
    k6 = k6_vec(i);
    Mx = [[-(k2+k3);k3;0],[0;-k5;k5],[k6;0;-k6]];
    Cx = concentration(k1,Mx,Ca,0,[0;0;0],t);
    Cxdata(:,i) = (alpha*Cx + Vb*Ca(t))';
end

figure('units','normalized','outerposition',[0 0 0.7 0.7]);
set(gcf,'Color','none');
hold on;
axis square;% grid on;
for i = 1:length(k6_vec)    
    if range(i) == 1
        plot(t,Cxdata(:,i), '-x', 'Markersize', 10, 'LineWidth', 3.5, ...
            'Displayname', strcat(['k6 = ', num2str(k6_vec(i))]));
    else
    plot(t,Cxdata(:,i), 'LineWidth', 3.5, 'Displayname', strcat(['k6 = ', num2str(k6_vec(i))]));   
    end
end
set(gca,'xtick',0:10:40); 
xlabel('time [min]','FontSize',30,'Interpreter','Latex'); 
ylabel('concentration $C_T$ [kBq/mL]','FontSize',30,'Interpreter','Latex');
legend('show', 'Location','northeastoutside')
ylim([0, 500])
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

fig = gcf; fig.PaperPositionMode = 'auto';
print(fullfile(path_figures, 'sensitivity_k6'),'-dpng','-r300')

%% Relative local sensitivity S_T
k1_or = 0.3;   % Values of the parameters provided by the numerical 
               % solution of the inverse problem.
k2_or = 0.4;
k3_or = 0.6;
k5_or = 0.6;
k6_or = 0.06;

t0 = 0;         % Inital data for the compartment model
C0 = [0; 0; 0];
y0 = [0; 0; 0];

n_par = 5;
range = [0.01, 0.1, 0.5, 1, 1.1, 1.5, 1.6, 1.7, 2, 10, 20, 100, 1000];
n_scaling_factor = numel(range);

S_T = zeros(n_par, n_scaling_factor);

%% k1
k2 = k2_or; k3 = k3_or; k5 = k5_or; k6 = k6_or;
for is = 1:n_scaling_factor
    k1 = k1_or*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t)', [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k1 = f_comp_sol_activity(A_coeff, Ca(t)', [1; 0; 0], t0, y0, t);
    ST_C_k1 = k1 * (y_k1 ./ C_time);
    
    CT_dk1 = alpha * y_k1;
    ST_k1 = abs(k1*(CT_dk1 ./ CT_time));
    
    S_T(1, is) = ST_k1(end);

end

clear k1 k2 k3 k5 k6

%% k2
k1 = k1_or; k3 = k3_or; k5 = k5_or; k6 = k6_or;
for is = 1:n_scaling_factor
    k2 = k2_or*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t)', [1; 0; 0], t0, C0, t);
    CT_time = (alpha * C_time + Vb * Ca(t) );
    
    y_k2 = f_comp_sol_activity(A_coeff, -C_time(1, :)', [1; 0; 0], t0, y0, t);
    ST_C_k2 = abs(k2 * (y_k2 ./ C_time));
    
    CT_dk2 = alpha * y_k2;
    ST_k2 = abs(k2 * CT_dk2 ./ CT_time);
    
    S_T(2, is) = ST_k2(end);

end

clear k1 k2 k3 k5 k6

%% k3 
k1 = k1_or; k2 = k2_or; k5 = k5_or; k6 = k6_or;
for is = 1:n_scaling_factor
    k3 = k3_or*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t)', [1; 0; 0], t0, C0, t);
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
k1 = k1_or; k2 = k2_or; k3 = k3_or; k6 = k6_or;
for is = 1:n_scaling_factor
    k5 = k5_or*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t)', [1; 0; 0], t0, C0, t);
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
k1 = k1_or; k2 = k2_or; k3 = k3_or; k5 = k5_or;
for is = 1:n_scaling_factor
    k6 = k6_or*range(is);

    A_coeff = [ [-(k2+k3); k3; 0], [0; -k5; k5], [k6; 0; -k6]];
    
    C_time = f_comp_sol_activity(A_coeff, k1*Ca(t)', [1; 0; 0], t0, C0, t);
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
hold on
plot(range, S_T(1, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^1$')
plot(range, S_T(2, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^2$')
plot(range, S_T(3, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^3$')
plot(range, S_T(4, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^5$')
plot(range, S_T(5, :), 'Linewidth', 3, 'Displayname',  '$\overline{S}_T^6$')
lgd = legend('show', 'Location', 'Bestoutside');
lgd.FontSize = 25;
set(lgd, 'Interpreter', 'Latex');
set(gca, 'XScale', 'log', 'FontSize', 25);
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ylabel('$\overline{S}_T^j$', 'FontSize', 25, 'Interpreter', 'Latex')
xlabel('Scaling factor', 'FontSize', 25, 'Interpreter', 'Latex')
axis tight
print(fullfile(path_figures, 'ST_all_parameters'), '-dpng','-r300')







