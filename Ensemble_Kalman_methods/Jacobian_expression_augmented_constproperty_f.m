% derive linearized Jacobian expression for augmented state vector [T, q]

% input: 
% state_augmented = [T_profile; q_magnitude]; % [T, q] is column vector of augmented state
% dt: time step for integration

% output: F is state transition matrix

function [F]...
    = Jacobian_expression_augmented_constproperty_f(state_augmented, dt)

%% geometry
ro_heater = 4.25*10^-3; % [m]
ro_zro2 = 30.70*10^-3; % [m] instead of 30.734 mm
dr = 0.05*10^-3; % [m] grid size: can be greater than greatest common divisor of two lengths

r_grid_heater = linspace(0, ro_heater, 1+ro_heater/dr)'; % including interface
r_grid_zro2 = linspace(ro_heater, ro_zro2, 1+(ro_zro2-ro_heater)/dr)'; % including interface
r_grid_total = linspace(0, ro_zro2, 1+ro_zro2/dr)';
N_grid_heater = length(r_grid_heater);
N_grid_total = length(r_grid_total);
I_heater = (r_grid_total < ro_heater) + 0; % 1 for heater, 0 for else

%% physical quantity
% dt = 1;
T_profile = state_augmented(1:end-1);
q_magnitude = state_augmented(end);

T_outer = T_profile(end-1); % [K] boundary condition

%% temperature dependent properties on grid
k_i = zeros(N_grid_total,1); % thermal conductivity at node i
k_i_face = zeros(N_grid_total-1,1); % correspond to i+1/2 location
rhocp_i = zeros(N_grid_total,1); % rho*c_p at node i

%% assign matrix values
F = zeros(N_grid_total + 1, N_grid_total + 1);
F_T = zeros(N_grid_total, N_grid_total);
F_q = zeros(N_grid_total + 1, 1);

%% thermal properties 
idx = 1:N_grid_heater-1; % heater
k_i(idx) = k_T_heater(T_profile(idx)); % use kelvin for input temperature
rhocp_i(idx) = rhocp_T_heater(T_profile(idx));

idx = N_grid_heater; % gap
k_i(idx) = 1/(0.5/k_T_heater(T_profile(idx)) + 0.5/k_T_zro2(T_profile(idx)) + Rgap_T(T_profile(idx))/dr);
rhocp_i(idx) = 0.5*(rhocp_T_heater(T_profile(idx)) + rhocp_T_zro2(T_profile(idx)));

idx = N_grid_heater+1:N_grid_total; % zro2
k_i(idx) = k_T_zro2(T_profile(idx));  % use kelvin for input temperature
rhocp_i(idx) = rhocp_T_zro2(T_profile(idx));

k_i_face = 1./(0.5./k_i(1:end-1) + 0.5./k_i(2:end));

%% F_T
% F_T symmetry
idx = 1;
F_T(idx,idx) = - 4*k_i_face(idx)/rhocp_i(idx);
F_T(idx,idx) = F_T(idx,idx)/dr/dr;

F_T(idx,idx+1) = 4*k_i_face(idx)/rhocp_i(idx);
F_T(idx,idx+1) = F_T(idx,idx+1)/dr/dr;

% F_T intermediate: heater, gap, insulation
idx = 2:(N_grid_total - 1);
vector_idx = idx + (idx - 1)*(N_grid_total);
F_T(vector_idx) = - (k_i_face(idx) + k_i_face(idx-1) + (k_i_face(idx) - k_i_face(idx-1))/2*dr./r_grid_total(idx) )./rhocp_i(idx);
F_T(vector_idx) = F_T(vector_idx)/dr/dr;

F_T(vector_idx+N_grid_total) = (1 + 1/2*dr./r_grid_total(idx)).*k_i_face(idx)./rhocp_i(idx);
F_T(vector_idx + N_grid_total) = F_T(vector_idx + N_grid_total)/dr/dr;

F_T(vector_idx- N_grid_total) = (1 - 1/2*dr./r_grid_total(idx)).*k_i_face(idx-1)./rhocp_i(idx);
F_T(vector_idx - N_grid_total) = F_T(vector_idx - N_grid_total)/dr/dr;

% F_T outer boundary
idx = N_grid_total;
F_T(idx, idx) = 0;

%% F_q
idx = 1:(N_grid_heater - 1);
F_q(idx) = 1./rhocp_i(idx);
F_q(N_grid_heater) = 1/2/rhocp_i(N_grid_heater);

%% calculation of F
F(1:end-1, 1:end-1) = F_T*dt + eye(N_grid_total);
F(1:end, end) = F_q*dt;
F(end, end) = 1; % q -> q

end
