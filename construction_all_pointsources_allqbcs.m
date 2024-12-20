% prepare operation for H (sensitivity) matrix for Gibbs sampling
% this code can be only run on desktop (memory issue)
clear;
close all;

% define input mass flux and heat source magnitude
G = 4215; % [kg/m2/s] 
heat_flux_source = 852*1000; % INPUT [W/m2]
T_inlet_uniform = 300; % [K] uniform temperature at inlet

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

%% physical properties
% water properties from XSteam
P_in = 1.00; % [bar] INPUT 'Pressure_inlet' from CDTS (12.00; 12.07; 12.04);
Tb_in = 26.2; % [degC]
k_cond = XSteam('tc_pT', P_in, Tb_in); % [W/m/K] thermal conductivity

%% load A matrix
digit_thousands = floor(G/1000);
digit_hundreds = floor((G - digit_thousands*1000)/100);
digit_tens= floor((G - digit_thousands*1000 - digit_hundreds*100)/10);
digit_ones= floor(G - digit_thousands*1000 - digit_hundreds*100 - digit_tens*10);
inputname      = ['./A_matrix_conv/A_matrix_conv_' num2str(digit_thousands) '' num2str(digit_hundreds) '' num2str(digit_tens) '' num2str(digit_ones) '.mat'];
load(inputname,'A_matrix');

%% matrix to save all temperatures and boundary conditions
T_steady_allqbcs = zeros(N_z*N_x_tot*N_y , (N_z-1)*N_y); % temperature (for all cases point heat source)
L_q1_allqbcs = sparse(N_z*N_x_tot*N_y, (N_z-1)*N_y); % thermal boundary condition (for all cases point heat source)

%% assign point heat source location
for source_idz=2:N_z
    for source_idy=1:N_y
        index_x = N_x_1;
        index_z = source_idz;
        index = index_z + (index_x - 1)*N_z + (source_idy - 1)*N_z*N_x_tot;
        
        id_qbcs = (source_idz-1) + (N_z-1)*(source_idy-1);
        L_q1_allqbcs(index, id_qbcs) = -(1)*dx*1/k_cond; % of point heat source of unit strength [W/m2]
    end
end

%% solve matrix equation for temperature (need large memory for computation)
T_steady_allqbcs = A_matrix\L_q1_allqbcs; % temperature map as column vector
T_steady_allqbcs = full(T_steady_allqbcs);

%% save result
% outputname = 'steady_3d_conduction_convection_turb_v1_G4215_Tqbcs.mat';
% save(outputname,'T_qbcs_allqbcs','-v7.3');
% fprintf('Result saved. \n')
