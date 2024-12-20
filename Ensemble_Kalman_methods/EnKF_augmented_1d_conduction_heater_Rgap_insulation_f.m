% Ensemble Kalman Filter implementation
% v2. make different time step (dt ~= 1)available
function [q_error_normalized, q_std_normalized, T_peak_error, T_peak_std]...
    = EnKF_augmented_1d_conduction_heater_Rgap_insulation_f(...
    sigma_noise, Rgap_noise, q_noise, dt, sensor_user_input, M_sample, inputname1)

%% load data: forward problem
% inputname1 = 'conduction1D_radial_Tdep_onelevel_qt_pulse1_Rgap_Tdep_dr005.mat';
load(inputname1,'time_grid','T_transient','R_transient','q_transient');

% original solutions
time_grid_solution = time_grid;
T_transient_solution = T_transient;
R_transient_solution = R_transient;
q_transient_solution = q_transient;

%% hyperparameter settings filtering method
rng('default');
% M_sample = 50; % INPUT
% sigma_noise = 0.1; % [K] INPUT noise magnitude
% q_noise = 10^8; % [W/m^3] INPUT 10^5
% Rgap_noise = 0.1; % [ratio] INPUT error physics modeling on Rgap

%% geometry
ro_heater = 4.25*10^-3; % [m]
ro_zro2 = 30.70*10^-3; % [m] instead of 30.734 mm
dr = 0.05*10^-3; % [m] grid size: can be greater than greatest common divisor of two lengths

r_grid_heater = linspace(0, ro_heater, 1+ro_heater/dr)'; % including interface
r_grid_zro2 = linspace(ro_heater, ro_zro2, 1+(ro_zro2-ro_heater)/dr)'; % including interface
r_grid_total = linspace(0, ro_zro2, 1+ro_zro2/dr)';
N_grid_heater = length(r_grid_heater);
N_grid_total = length(r_grid_total);

%% temperature sensor placement (index)
% sensor_user_input = [146];

Nsensor = length(sensor_user_input);
Hobs = zeros(Nsensor, N_grid_total + 1); % observation matrix for augmented state [T, q]
for i = 1:Nsensor % Fill in the observation matrix
    Hobs(i, sensor_user_input(i)) = 1;
end
Gamma = diag(sigma_noise * ones(1, Nsensor)); % noise covariance matrix

%% Inverse time steps
time_end = 100; % [s] INPUT last time
% dt = 2; % [s] INPUT inverse time step size
time_grid = linspace(0, time_end, 1+time_end/dt)';
N_time = length(time_grid);

%% translate original data onto inverse time step
T_transient_solution = myLinearInterp(time_grid_solution, T_transient', time_grid)';
R_transient_solution = myLinearInterp(time_grid_solution, R_transient, time_grid);
q_transient_solution = myLinearInterp(time_grid_solution, q_transient, time_grid);
noise = sigma_noise * randn([size(T_transient_solution) M_sample]); % Generating Gaussian noise
T_transient_solution_noisy =  repmat(T_transient_solution, [1, 1, M_sample]) + noise; % Adding the noise to the data matrix

%% matrix (for individual sample): state, thermal properties, data
T_transient_sample = zeros(N_grid_total, N_time, M_sample);
q_transient_sample = zeros(N_time, M_sample);
R_transient_sample = zeros(N_time, M_sample); % [m^2-K/W]
Yobs_sample = zeros(Nsensor, N_time, M_sample);
Yobs_sample = repmat(Hobs(:,1:end-1)*T_transient_solution_noisy(:,1:N_time), [1, 1, M_sample]) + sigma_noise*randn(Nsensor, N_time, M_sample); % corrupted measurement

k_i = zeros(N_grid_total, M_sample); % thermal conductivity at node i
rhocp_i = zeros(N_grid_total, M_sample); % rho*c_p at node i
k_i_face = zeros(N_grid_total-1, M_sample); % correspond to i+1/2 location

%% mean and std (estimation = after Kamlan correction)
T_transient_est_mean = zeros(N_grid_total, N_time);
q_transient_est_mean = zeros(N_time, 1);
R_transient_est_mean = zeros(N_time, 1);
T_transient_est_std = zeros(N_grid_total, N_time);
q_transient_est_std = zeros(N_time, 1);
R_transient_est_std = zeros(N_time, 1);

%% statistics  
Cov_MC_temp = zeros(N_grid_total+1, N_grid_total+1); % estimated covariance (updated each time step)
Var_MC = zeros(N_grid_total+1, N_time); % store variance of each state
K_MC = zeros(N_grid_total+1, Nsensor, N_time); % Kalman gain
state_augmented_sample_temp = zeros(N_grid_total+1, M_sample); 

%% physical properties
T_outer = 500 + 273.15; % [K] boundary condition

%% initialization (t=0)
it = 1;
% fprintf('initialization (time: %d second) \n', time_grid(it));

T_transient_sample(:,it,:) = T_outer;
state_augmented_sample_temp(1:end-1, :) = T_transient_sample(:,it,:);
state_augmented_sample_temp(end, :) = q_transient_sample(it,:);
Cov_MC_temp = cov(state_augmented_sample_temp');

idx = 1:N_grid_heater-1; % heater
k_i(idx,:) = k_T_heater(T_outer); % use kelvin for input temperature
rhocp_i(idx,:) = rhocp_T_heater(T_outer);

idx = N_grid_heater; % gap
R_transient_sample(it,:) = Rgap_T(T_outer)*(1 + Rgap_noise*randn(1,M_sample));
k_i(idx,:) = 1./(0.5/k_T_heater(T_outer) + 0.5/k_T_zro2(T_outer) + R_transient_sample(it,:)/dr);
rhocp_i(idx,:) = 0.5*(rhocp_T_heater(T_outer) + rhocp_T_zro2(T_outer));

idx = N_grid_heater+1:N_grid_total; % zro2
k_i(idx,:) = k_T_zro2(T_outer);
rhocp_i(idx,:) = rhocp_T_zro2(T_outer);

k_i_face = 1./(0.5./k_i(1:end-1,:) + 0.5./k_i(2:end, :));
% face #1: (1,2), face #2: (2,3), face #(N_grid_total-1):(N_grid_total-1, N_grid_total)

Var_MC(:, it) = sqrt(diag(Cov_MC_temp));

T_transient_est_mean(:,it) =  mean(T_transient_sample(:,it,:), 3);
q_transient_est_mean(it) = mean(q_transient_sample(it,:));
R_transient_est_mean(it) = mean(R_transient_sample(it,:));
T_transient_est_std(:,it) =  std(T_transient_sample(:,it,:), 0, 3);
q_transient_est_std(it) = std(q_transient_sample(it,:), 0);
R_transient_est_std(it) = std(R_transient_sample(it,:), 0);

%% solve inverse problem (t>0)
for it=2:N_time
%     fprintf('time: %d second \n', time_grid(it));
    
    %% evolution heat source (process noise)
    q_transient_sample(it,:) = q_transient_sample(it-1, :) + q_noise*randn(1,M_sample);
    
for m=1:M_sample % for each sample

    %% assign matrix values
    A = sparse(N_grid_total, N_grid_total); % new time step
    B = sparse(N_grid_total, N_grid_total); % prev time step
    C = sparse(N_grid_total, 1); % source term
    
    %% discretization for T equation
    % center symmetry line
    idx = 1;
    A(idx, idx) = 4*k_i_face(idx,m) + rhocp_i(idx,m)/dt*dr^2;
    A(idx, idx+1) = -4*k_i_face(idx,m);
    B(idx, idx) = rhocp_i(idx,m)/dt*dr^2;
    C(idx) = q_transient_sample(it, m)*dr^2; % keep using the one from previous time step

    % intermediate
    idx = 2:(N_grid_total - 1);
    vector_idx = idx + (idx - 1)*N_grid_total;
    A(vector_idx) = (k_i_face(idx,m)+k_i_face(idx-1,m)) + (k_i_face(idx,m)-k_i_face(idx-1,m))/2./r_grid_total(idx)*dr + rhocp_i(idx,m)/dt*dr^2; % [i, i]
    A(vector_idx+N_grid_total) = -k_i_face(idx,m) - k_i_face(idx,m)/2./r_grid_total(idx)*dr; % [i, i+1]
    A(vector_idx-N_grid_total) = -k_i_face(idx-1,m) + k_i_face(idx-1,m)/2./r_grid_total(idx)*dr; % [i, i-1]
    B(vector_idx) = rhocp_i(idx,m)/dt*dr^2;
    
    idx = 2:N_grid_heater-1; % heat source from heater (excluding gap)
    C(idx) = q_transient_sample(it, m)*dr^2;
    C(N_grid_heater) = 1/2*q_transient_sample(it, m)*dr^2;
    
    % outer boundary
    idx = N_grid_total;
    A(idx, idx) = 1;
    C(idx) = T_outer;    
    
    % solve transient equation
    T_transient_sample(:,it,m) = A\(B*T_transient_sample(:,it-1,m) + C);
    
end

% augmented state
state_augmented_sample_temp(1:end-1, :) = T_transient_sample(:,it,:);
state_augmented_sample_temp(end, :) = q_transient_sample(it,:);
Cov_MC_temp = cov(state_augmented_sample_temp');

% Kalman gain and update
K_MC(:,:,it) = Cov_MC_temp*Hobs'*inv(Gamma + Hobs*Cov_MC_temp*Hobs');

% assimilate individual result
for m=1:M_sample % for each sample
state_augmented_sample_temp(:,m) = ...
    state_augmented_sample_temp(:,m) + K_MC(:,:,it)*(Yobs_sample(:,it,m) - Hobs*state_augmented_sample_temp(:,m));
end

% assign values
T_transient_sample(:,it,:) = state_augmented_sample_temp(1:end-1, :);
q_transient_sample(it, :) = state_augmented_sample_temp(end, :);

%% properties for next time step (one-level scheme)
for m=1:M_sample % for each sample    
    idx = 1:N_grid_heater-1; % heater
    k_i(idx,m) = k_T_heater(T_transient_sample(idx,it,m)); % use kelvin for input temperature
    rhocp_i(idx,m) = rhocp_T_heater(T_transient_sample(idx,it,m));

    idx = N_grid_heater; % gap
    R_transient_sample(it,m) = Rgap_T(T_transient_sample(idx,it,m))*(1 + Rgap_noise*randn());
    k_i(idx,m) = 1/(0.5/k_T_heater(T_transient_sample(idx,it,m)) + 0.5/k_T_zro2(T_transient_sample(idx,it,m)) + R_transient_sample(it,m)/dr);
    rhocp_i(idx,m) = 0.5*( rhocp_T_heater(T_transient_sample(idx,it,m)) + rhocp_T_zro2(T_transient_sample(idx,it,m)) );

    idx = N_grid_heater+1:N_grid_total; % zro2
    k_i(idx,m) = k_T_zro2(T_transient_sample(idx,it,m));  % use kelvin for input temperature
    rhocp_i(idx,m) = rhocp_T_zro2(T_transient_sample(idx,it,m));

    k_i_face(:,m) = 1./(0.5./k_i(1:end-1,m) + 0.5./k_i(2:end,m));
end

% mean and std 
T_transient_est_mean(:,it) = mean(T_transient_sample(:,it,:), 3);
q_transient_est_mean(it) = mean(q_transient_sample(it,:));
R_transient_est_mean(it) = mean(R_transient_sample(it,:));
T_transient_est_std(:,it) = std(T_transient_sample(:,it,:), 0, 3);
q_transient_est_std(it) = std(q_transient_sample(it,:), 0);
R_transient_est_std(it) = std(R_transient_sample(it,:), 0);

end

%% calculate reconstruction errors
id_plot = 1:N_time-1; % for all inverse time steps

% heat source
q_transient_diff = q_transient_est_mean(id_plot) - q_transient_solution(id_plot); % discrepancy vector
q_error_normalized = mean(abs(q_transient_diff))/mean(q_transient_solution); % bias (normalized)
q_std_normalized = mean(q_transient_est_std)/mean(q_transient_solution); % std (normalized)

% peak temperature
id_peak = 1;
T_peak_error = mean(abs(T_transient_est_mean(id_peak,id_plot) - T_transient_solution(id_peak,id_plot))); % bias
T_peak_std = mean(T_transient_est_std(id_peak,id_plot)); % std

% fprintf('Δq/q0 = %.2f ± %.2f \n', q_error_normalized, q_std_normalized);
% fprintf('ΔT(peak) = %.2f ± %.0f [K] \n', T_peak_error, T_peak_std);

end