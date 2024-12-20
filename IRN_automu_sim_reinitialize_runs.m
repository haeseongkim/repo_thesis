% to plot convergence (statistics) result with different sensor placements
clear;
close all;

mu0 = 0.1405;

noise_magnitude_case = [0.1, 0.5, 1.0]; % INPUT
N_noise_case = length(noise_magnitude_case);

input_z_case = {}; % INPUT below
input_y_case = {}; % INPUT below
input_z_case{1} = [2:1:101];
input_y_case{1} = [1:1:19];
input_z_case{2} = [2:2:100];
input_y_case{2} = [1:2:19];
input_z_case{3} = [2:3:101];
input_y_case{3} = [1:3:19];
input_z_case{4} = [2:4:98];
input_y_case{4} = [1:4:17];
input_z_case{5} = [2:5:97];
input_y_case{5} = [1:5:16];
input_z_case{6} = [2:6:98];
input_y_case{6} = [1:6:19];

sensor_distances = [1, 2, 3, 4, 5, 6];
N_sensor_case = length(input_z_case);

mean_qbcs_reconst_diff_save = zeros(N_noise_case, N_sensor_case); % same row: same noise magnitude
reconst_Lqq_save = zeros(N_noise_case, N_sensor_case); 
reconst_nnz_save = zeros(N_noise_case, N_sensor_case);
noise_estimation_save = zeros(N_noise_case, N_sensor_case);
iou_save = zeros(N_noise_case, N_sensor_case);
solution_Lqq = 0;
solution_L0 = 0;

%% input data
% inputname1 = './single_qsource/steady_3d_conduction_convection_turb_v1_single_qsource_rectangle_G4215_0852kW.mat';
inputname1 = './single_qsource/steady_3d_conduction_convection_turb_v1_single_qsource_halfcross_G4215_0852kW.mat';
% inputname1 = './single_qsource/steady_3d_conduction_convection_turb_v1_single_qsource_arc1_G4215_0852kW.mat';
% inputname1 = './single_qsource/steady_3d_conduction_convection_turb_v1_single_qsource_groove1_G4215_0852kW.mat';

%% calculation of each case
% comment the calculation loop out when we have data
for idn=1:N_noise_case
    for ids=1:N_sensor_case
        [mean_qbcs_reconst_diff_save(idn,ids), reconst_Lqq_save(idn,ids), reconst_nnz_save(idn,ids), ...
            noise_estimation_save(idn,ids), iou_save(idn,ids), solution_Lqq, solution_L0]...
        = IRN_automu_sim_reinitialize_f(noise_magnitude_case(idn), input_z_case{ids}, input_y_case{ids}, mu0, inputname1);    
    end
end

%% load pre-calculated result if available 
% load('IRN_automu_runsf2_rectangle.mat');
% load('IRN_automu_runsf2_halfcross.mat');
% load('IRN_automu_runsf2_arc1.mat');
% load('IRN_automu_runsf2_groove1.mat');

%% plot result
fontsize = 20;
figure(1); 
set(gcf, 'WindowState', 'maximized');
%% MAE
subplot(3,2,1)
for idn=1:N_noise_case
    plot(sensor_distances, mean_qbcs_reconst_diff_save(idn,:),'-o','Linewidth',2,'displayname', sprintf(' noise %.1f K', noise_magnitude_case(idn)))
    hold on
end

xlabel('sensor distance [pixel]');
ylabel('Heat flux error [W/m^2]')
title('Mean absolute error |x_{recon} - x_{sol}|')
set(gca,'FontSize',fontsize,'fontname','times')
set(gca, 'YScale', 'log')
xlim([1 6])
ylim([10^1 10^5])
yticks([10^1 10^2, 10^3, 10^4, 10^5]);
% yticklabels({'10^{-3}', '10^{-2}', '10^{-1}', '1', '10'});
grid on; 
grid minor;

%% Regularization
subplot(3,2,3)
for idn=1:N_noise_case
    plot(sensor_distances, reconst_Lqq_save(idn,:),'-o','Linewidth',2,'displayname', sprintf(' noise %.1f K', noise_magnitude_case(idn)))
    hold on
end
yline(solution_Lqq, '--','LineWidth', 2,'displayname','solution')

xlabel('sensor distance [pixel]');
ylabel('$\|D\mathbf{x}\|_q^q$', 'Interpreter', 'latex', 'FontSize', fontsize);
legend('show','location','southeast')
title('Regularization');
set(gca,'FontSize',fontsize,'fontname','times')
set(gca, 'YScale', 'log')
ylim([10 3*10^3])
yticks([10^1, 10^2, 10^3 3*10^3]);
grid on; 
grid minor;

%% number edge
subplot(3,2,2)
for idn=1:N_noise_case
    plot(sensor_distances, reconst_nnz_save(idn,:),'-o','Linewidth',2,'displayname', sprintf(' noise %.1f K', noise_magnitude_case(idn)))
    hold on
end
yline(solution_L0,'--','LineWidth',2,'displayname',' solution')

xlabel('sensor distance [pixel]');
ylabel('$\|Lx\|_0$', 'Interpreter', 'latex', 'FontSize', fontsize);
title('Edge feature');
set(gca,'FontSize', fontsize,'fontname','times')
set(gca, 'YScale', 'log')
ylim([50 3*10^3])
yticks([10^2, 10^3 3*10^3]);
grid on; 
grid minor;

%% intersection over union
subplot(3,2,4)
for idn=1:N_noise_case
    plot(sensor_distances, iou_save(idn,:),'-o','Linewidth',2,'displayname', sprintf(' noise %.1f K', noise_magnitude_case(idn)))
    hold on
end
yline(1,'--','Linewidth',2,'handlevisibility','off')

xlabel('sensor distance [pixel]');
ylabel('Intersection over Union');
ylim([0 1.2])
title('Area feature')
set(gca, 'FontSize', fontsize, 'FontName', 'Times');
grid on
grid minor

%% estimated noise
subplot(3,2,5)
for idn=1:N_noise_case
    plot(sensor_distances, noise_estimation_save(idn,:),'-o','Linewidth',2,'displayname', sprintf(' noise %.1f K', noise_magnitude_case(idn)))
    hold on
    yline(noise_magnitude_case(idn),'--','Linewidth',2,'handlevisibility','off')
    hold on
end

xlabel('sensor distance [pixel]');
ylabel('Noise magnitude [K]');
title('Estimated noise')
set(gca, 'FontSize', fontsize, 'FontName', 'Times');
grid on
grid minor

%% plot original solution
load(inputname1,'z_grid_mm','y_grid_mm','q_source','vertices');

%% produce sensor matrix 
input_z = [2:1:101];
input_y = [1:1:19];
sensor_user_input = zeros(length(input_z)*length(input_y),3);
for iz=1:length(input_z)
    for iy=1:length(input_y)
        row = (iz-1)*length(input_y) + iy;
        sensor_user_input(row, 1) = input_z(iz);
        sensor_user_input(row, 3) = input_y(iy);
    end
end
sensor_index_z = sensor_user_input(:,1);
sensor_index_y = sensor_user_input(:,3);

figure(1)
subplot(3,2,6)
imagesc(z_grid_mm(2:end), y_grid_mm, q_source(2:end,:)')
hold on
scatter(z_grid_mm(sensor_index_z), y_grid_mm(sensor_index_y),4,'displayname','sensor locations','MarkerFaceColor','#EDB120')
hold on
plot(vertices(:,1)*1000, vertices(:,2)*1000, 'w--', 'LineWidth',2,'handlevisibility','off'); % draw white line (need existing figure)
colormap('jet');
c = colorbar; % Get the handle of the colorbar
ylabel(c, 'heat flux W/m^2')
xlabel('z grid [mm]')
ylabel('y grid [mm]')
titleStr = sprintf('original solution');
title(titleStr);
legend('show','location','southeast')
set(gca,'FontSize',fontsize,'fontname','times')
grid on
grid minor

