% pre-computation of division by A_matrix (takes > 5 minutes for A\b)
clear;

G = 4215; % INPUT
T_inlet_uniform = 300; % [K]

%% load A matrix
digit_thousands = floor(G/1000);
digit_hundreds = floor((G - digit_thousands*1000)/100);
digit_tens= floor((G - digit_thousands*1000 - digit_hundreds*100)/10);
digit_ones= floor(G - digit_thousands*1000 - digit_hundreds*100 - digit_tens*10);
inputname      = ['./A_matrix_conv/A_matrix_conv_' num2str(digit_thousands) '' num2str(digit_hundreds) '' num2str(digit_tens) '' num2str(digit_ones) '.mat'];
load(inputname,'A_matrix','A_conv');

%% number grid
e_channel = 0.5*1.95*10^-3; % [m] half height of flow channel (we use symmetry here: total channel height is 2 mm)
L_channel = 0.056; % [m] flow length
w_channel = 18*10^-3; % [m] width channel (just for information)
% T_inlet_uniform = T_K(1); % [K] uniform temperature at inlet
t_solid = 2.00*10^-3; % [m] solid thickness

t_liquid_boundary = (0.4875*10^-6)*(50); % [m] boundary layer thickness: (near-wall)*(# cells)
N_x_1 = 300+1; % INPUT total number of liquid grid points (including both ends) in x-direction
N_x_1_boundary = 50+1; % INPUT including liquid-solid interface
N_x_1_bulk = N_x_1 - N_x_1_boundary; % including bulk end, excluding boundary
N_z = 100+1; % INPUT number of grid points (including both ends: inlet and outlet) in z-direction
N_x_2 = 40+1; % INPUT for solid substrate (including both ends: interface and outer wall)
N_x_tot = N_x_1 + N_x_2 - 1;
N_y = 19; % INPUT number of grid points in width direction

% First point is under specified temperature, and the rest (N_z - 1) points will be under heat flux.
z_grid = linspace(0, L_channel, N_z); % [m] distance from inlet
z_grid_mm = z_grid*1000;
dz = L_channel/(N_z - 1); % equal grid size
y_grid = linspace(0, w_channel/2, N_y); % [m] distance from the symmetry plane
y_grid_mm = y_grid*1000;
dy = y_grid(2) - y_grid(1); 
x_grid_liquid_boundary = linspace(e_channel-t_liquid_boundary, e_channel, N_x_1_boundary);

% hyperbolic tangent grid with stretching parameter
beta_stretch = 2.1; % INPUT, larger the steeper at the end. check smoothness of the cell size
% Create the one-sided hyperbolic-tangent spaced grid
xi = linspace(0, 1, N_x_1_bulk + 1);
grid_hyperbolic = x_grid_liquid_boundary(1) - (1-tanh(beta_stretch*xi)/tanh(beta_stretch))*(e_channel - t_liquid_boundary); % including overlap
x_grid_liquid_bulk = grid_hyperbolic(1:end-1); % excluding interface, sees non-uniform grid

x_grid_liquid = [x_grid_liquid_bulk, x_grid_liquid_boundary]; % [m] distance from bulk (row vector)
dx = x_grid_liquid(end)-x_grid_liquid(end-1); % near-wall cell size
x_grid_solid = linspace(e_channel+dx, e_channel+t_solid, N_x_2-1); % excluding interface
dx_solid = (x_grid_solid(end)-x_grid_solid(1))/(N_x_2 - 2);

x_grid_total = [x_grid_liquid, x_grid_solid]; % entire grid

%% input heat flux shape: constant
q_source1 = zeros(N_z, N_y);
heater_start_z = (55.95/2-19.5/2)*10^-3; % INPUT [m]
heater_end_z = (55.95/2+19.5/2)*10^-3; % INPUT [m]
heater_start_y = 0*10^-3; % INPUT [m]
heater_end_y = 4.5*10^-3 + sqrt(eps); % INPUT [m]

q_location_z = (z_grid <= heater_end_z & z_grid >= heater_start_z); % around middle
q_location_y = (y_grid <= heater_end_y & y_grid >= heater_start_y); 
q_source1 = 1*q_location_z'*q_location_y; % 2d distribution (heat_flux_source = 1 W/m2)

%% properties
P_in = 1.00; % [bar] INPUT 'Pressure_inlet' from CDTS (12.00; 12.07; 12.04);
Tb_in = 26.2; % [degC]
k_cond = XSteam('tc_pT', P_in, Tb_in); % [W/m/K] thermal conductivity

%% heat source 
L_qbcs1 = zeros(N_z*N_x_tot*N_y, 1); % for unit heat source
for index_y=1:N_y
    % interface: index_x = N_x_1, index_z = N_z (node d)
    index_x = N_x_1;
    index_z = N_z;
    index = index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;
    L_qbcs1(index) = -(1)*dx*q_source1(index_z, index_y)/k_cond;

    % interface: index_x = N_x_1 (node c)
    vector_index_z = 2:(N_z-1);
    index_x = N_x_1;
    vector_index = vector_index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;
    L_qbcs1(vector_index) = -(1)*dx*q_source1(vector_index_z + (index_y-1)*N_z)/k_cond;   
end

%% inlet temperature
L_Tbcs = sparse(N_z*N_x_tot*N_y,1);
for index_y=1:N_y
    % inlet: index_z = 1 for liquid and solid
    index_z = 1;
    vector_index_x = 1:N_x_tot;
    vector_index = index_z + (vector_index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;
    L_Tbcs(vector_index) = T_inlet_uniform;
end

%% division with A matrix
L_merged = [L_qbcs1, L_Tbcs];

% div1
A_div_L_merged = A_matrix\L_merged;
% div2
A_div_A_conv_A_div_L_merged = A_matrix\(A_conv*A_div_L_merged);

%% Results for save: all grid
A_div_L_qbcs1 = A_div_L_merged(:,1);
A_div_A_conv_A_div_L_qbcs1 = A_div_A_conv_A_div_L_merged(:,1);
A_div_A_conv_A_div_L_Tbcs = A_div_A_conv_A_div_L_merged(:,2);

%% sensor location (use index_y=1): all sensors centerline (including inlet)
sensor_user_input = zeros(N_z,3); sensor_user_input(:,1) = 1:N_z; sensor_user_input(:,2) = N_x_1; sensor_user_input(:,3) = 1; % sensors all grid
% produce sensor matrix
sensor_index_z = sensor_user_input(:,1);
sensor_index_x = sensor_user_input(:,2);
sensor_index_y = sensor_user_input(:,3);
number_sensor = length(sensor_index_z);
sensor_index = sensor_index_z + (sensor_index_x-1)*N_z + (sensor_index_y-1)*N_z*N_x_tot; %

Phi_location_sensor = sparse(number_sensor, N_x_tot*N_z*N_y);
for ids=1:size(sensor_index, 1)    
    Phi_location_sensor(ids, sensor_index(ids)) = 1;
end

%% Results for save: centerline data
centerline_A_div_L_qbcs1 = Phi_location_sensor*A_div_L_qbcs1;
centerline_A_div_L_qbcs1 = full(centerline_A_div_L_qbcs1);

centerline_A_div_A_conv_A_div_L_qbcs1 = Phi_location_sensor*A_div_A_conv_A_div_L_qbcs1;
centerline_A_div_A_conv_A_div_L_qbcs1 = full(centerline_A_div_A_conv_A_div_L_qbcs1);

centerline_A_div_A_conv_A_div_L_Tbcs = Phi_location_sensor*A_div_A_conv_A_div_L_Tbcs;
centerline_A_div_A_conv_A_div_L_Tbcs = full(centerline_A_div_A_conv_A_div_L_Tbcs);
