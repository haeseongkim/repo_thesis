function [A] = construction_Amatrix_v1(G)
% input: G (mass flow rate), output: discretization matrix A

e_channel = 0.5*1.95*10^-3; % [m] half height of flow channel (we use symmetry here: total channel height is 2 mm)
L_channel = 0.056; % [m] flow length
w_channel = 18*10^-3; % [m] width channel (just for information)
T_inlet_uniform = 300;
t_solid = 2.00*10^-3; % [m] solid thickness

t_liquid_boundary = (0.4875*10^-6)*(50); % INPUT [m] boundary layer thickness: (near-wall)*(# cells)
N_x_1 = 300+1; % INPUT total number of liquid grid points (including both ends) in x-direction
N_x_1_boundary = 50+1; % INPUT including liquid-solid interface
N_x_1_bulk = N_x_1 - N_x_1_boundary; % including bulk end, excluding boundary
N_z = 100+1; % INPUT number of grid points (including both ends: inlet and outlet) in z-direction
N_x_2 = 40+1; % INPUT for solid substrate (including both ends: interface and outer wall)
N_x_tot = N_x_1 + N_x_2 - 1;
N_y = 10; % INPUT number of grid points in width direction

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
ratio_xz = dx/dz;
x_grid_solid = linspace(e_channel+dx, e_channel+t_solid, N_x_2-1); % excluding interface
dx_solid = (x_grid_solid(end)-x_grid_solid(1))/(N_x_2 - 2);
ratio_xs = dx/dx_solid;

% % uniform grid
x_grid_solid = linspace(e_channel+dx, e_channel+t_solid, N_x_2-1);

x_grid_total = [x_grid_liquid, x_grid_solid]; % entire grid

%% input heat flux shape: constant
q_source = zeros(N_z, N_y);
heater_start_z = (55.95/2-19.5/2)*10^-3; % INPUT [m]
heater_end_z = (55.95/2+19.5/2)*10^-3; % INPUT [m]
heater_start_y = 0*10^-3; % INPUT [m]
heater_end_y = 4.5*10^-3 + sqrt(eps); % INPUT [m]
heat_flux_source = 852*1000; % INPUT [W/m2] 585(02), 848(03), 852(04), 1546(07), 2224(10)

q_location_z = (z_grid <= heater_end_z & z_grid >= heater_start_z); % around middle
q_location_y = (y_grid <= heater_end_y & y_grid >= heater_start_y); 
q_source = heat_flux_source*q_location_z'*q_location_y; % 2d distribution

%% physical properties
% water properties from XSteam
P_in = 1.00; % [bar] INPUT 'Pressure_inlet' from CDTS (12.00; 12.07; 12.04);
Tb_in = 26.2; % [degC]

rho = XSteam('rho_pT', P_in, Tb_in); % [kg/m3] density
mu = XSteam('my_pT', P_in, Tb_in); % [Pa-s] dynamic viscosity
k_cond = XSteam('tc_pT', P_in, Tb_in); % [W/m/K] thermal conductivity
cp = XSteam('Cp_pT', P_in, Tb_in)*1000; % [J/kg] heat capacity
Pr = mu*cp/k_cond;
nu = mu/rho;
alpha = k_cond/(rho*cp);

% sapphire properties
rho_s = 3980.0; % [kg/m^3]
k_cond_s = 35; % [W/m/K]
cp_s = 929.0; % [J/kg-K]
alpha_s = k_cond_s/(rho_s*cp_s);

% interface properties
k_cond_int = 2*k_cond_s*k_cond/(k_cond_s+k_cond); % at interface
k_cond_lint = 4*k_cond_s*k_cond/(3*k_cond_s+1*k_cond); % between liquid and interface
k_cond_sint = 4*k_cond_s*k_cond/(1*k_cond_s+3*k_cond); % between solid and interface

%% dimensionless numbers
% velocity profile expression for fully developed turbulent pipe flow
vz_avg = G/rho; % [m/s] average flow velocity
D_h = 4*(2*e_channel*w_channel)/(4*e_channel+2*w_channel); % [m] hydraulic diameter of 3D channel, use it when needed
Re_b = G*D_h/mu;
Nu = 0.023*Re_b^0.8*Pr^0.4; % Dittus-Boelter correlation, heating
htc_DB = Nu*k_cond/(2*D_h);
tau_wall = 0.0332*rho*vz_avg*vz_avg*(mu/(rho*D_h*0.5*vz_avg))^(1/4); % wall shear stress, derived from Blasius empirical formulation (Re 4000~10^5)
v_friction = sqrt(tau_wall/rho); % [m/s] friction velocity

Pr_turb = 0.9; % [-] turbulent Prandtl number = viscous/thermal, 0.9 is the simplest approach.
vt = 1.15; % [-] multiplication for velocity magnitude (2D velocity field for 3D simulation)

%%
vz_profile = zeros(N_x_1,1); % axial velocity profile, function of x
eps_H = zeros(N_x_1,1); % [m2/s] eddy diffusivity of heat
T_inlet_developed = zeros(N_x_tot,1); % temperature profile at inlet

%%  inlet profiles: velocity, inlet temperature, eddy diffusivity of heat
for index_x = 1:N_x_tot
    T_inlet_developed(index_x) = T_inlet_uniform; % uniform inlet temperature
end
for index_x = 1:N_x_1
        vz_profile(index_x) = vt*velocity_profile_turb(x_grid_liquid(index_x), e_channel, v_friction, nu); % x is distance from bulk so makes sense
        eps_H(index_x) = (1/Pr_turb)*eddy_diffusivity_momentun_VanDriest(x_grid_liquid(index_x), e_channel, v_friction, nu); % eddy diffusivity of heat
end

%% A matrix
A = sparse(N_z*N_x_tot*N_y, N_z*N_x_tot*N_y);
% b = zeros(N_z*N_x_tot*N_y, 1);
% I_location_qbcs = sparse(N_z*N_x_tot*N_y, N_z*N_x_tot*N_y); % imposed heat flux boundary condition
% I_location_Tbcs = sparse(N_z*N_x_tot*N_y, N_z*N_x_tot*N_y); % display locations of temperature boundary condition

%% assign values on matrix. 
% note: (-1)*(1/k_liquid)*dx^2 is multiplied in equation to avoid division by small number.
for index_y=1:N_y
fprintf('index_y = %d (N_y=%d) \n', index_y, N_y);
% liquid
%% inlet: index_z = 1
index_z = 1;
vector_index_x = 1:N_x_1;
vector_index = index_z + (vector_index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = 1; % (i, j, k)
% b(vector_index) = T_inlet_developed(vector_index_x);
% I_location_Tbcs(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = 1; % (i, j, k)

%% outermost wall: index_x = 1, index_z = N_z (node e)
index_z = N_z;
index_x = 1;
index = index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_e = x_grid_liquid(index_x+1) - x_grid_liquid(index_x); % bulk end
dx_w = dx_e; % symmetry
k_turb_P = rho*cp*eps_H(index_x);
k_turb_E = rho*cp*eps_H(index_x+1);
k_turb_W = k_turb_E; % symmetry
k_turb_e = 2*k_turb_P*k_turb_E/(k_turb_P+k_turb_E);
k_turb_w = 2*k_turb_P*k_turb_W/(k_turb_P+k_turb_W);
k_total_P = k_cond + k_turb_P; % constant k_cond
k_total_e = k_cond + k_turb_e; 
k_total_w = k_cond + k_turb_w; % symmetry
v_P = vz_profile(index_x);

a_w = 2*(k_total_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_total_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_total_P/dz/dz; % north, i+1
a_s = k_total_P/dz/dz + rho*cp*v_P/dz; % south, i-1
a_i = k_total_P/dy/dy; % inward, k-1
a_o = k_total_P/dy/dy; % outward, k+1

A(index, index) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
% (i, j-1, k)
A(index, index + N_z) = -dx*dx/k_cond*(-a_w-a_e); % (i, j+1, k)
A(index, index - 1) = -dx*dx/k_cond*(-a_n - a_s); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_o); % (i, j, k+1)
end

%% outermost wall: index_x = 1 (node f)
vector_index_z = 2:(N_z-1);
index_x = 1;
vector_index = vector_index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_e = x_grid_liquid(index_x+1) - x_grid_liquid(index_x); % bulk end
dx_w = dx_e; % symmetry
k_turb_P = rho*cp*eps_H(index_x);
k_turb_E = rho*cp*eps_H(index_x+1);
k_turb_W = k_turb_E; % symmetry
k_turb_e = 2*k_turb_P*k_turb_E/(k_turb_P+k_turb_E);
k_turb_w = 2*k_turb_P*k_turb_W/(k_turb_P+k_turb_W);
k_total_P = k_cond + k_turb_P; % constant k_cond
k_total_e = k_cond + k_turb_e; 
k_total_w = k_cond + k_turb_w; % symmetry
v_P = vz_profile(index_x);

a_w = 2*(k_total_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_total_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_total_P/dz/dz; % north, i+1
a_s = k_total_P/dz/dz + rho*cp*v_P/dz; % south, i-1
a_i = k_total_P/dy/dy; % inward, k-1
a_o = k_total_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i. j, k)
% (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_w-a_e);% (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

%% near interface: index_x=N_x_1 - 1, index_z = N_z (node bd)
index_z = N_z;
index_x = N_x_1-1;
index = index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_e = x_grid_liquid(index_x+1) - x_grid_liquid(index_x); % bulk end
dx_w = x_grid_liquid(index_x) - x_grid_liquid(index_x-1); % 
k_turb_P = rho*cp*eps_H(index_x);
k_turb_E = 0; % rho*cp*eps_H(index_x+1)
k_turb_W = rho*cp*eps_H(index_x-1); 
k_turb_e = 0; % 2*k_turb_P*k_turb_E/(k_turb_P+k_turb_E);
k_turb_w = 2*k_turb_P*k_turb_W/(k_turb_P+k_turb_W);
k_cond_e = 4*k_cond_s*k_cond/(3*k_cond_s+1*k_cond); % between liquid and interface
k_cond_P = k_cond;
k_cond_w = k_cond;
k_total_P = k_cond_P + k_turb_P; % 
k_total_e = k_cond_e + k_turb_e; 
k_total_w = k_cond_w + k_turb_w; %
v_P = vz_profile(index_x);

a_w = 2*(k_total_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_total_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_total_P/dz/dz; % north, i+1
a_s = k_total_P/dz/dz + rho*cp*v_P/dz; % south, i-1
a_i = k_total_P/dy/dy; % inward, k-1
a_o = k_total_P/dy/dy; % outward, k+1

A(index, index) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(index, index - N_z) = -dx*dx/k_cond*(-a_w); % (i, j-1, k)
A(index, index + N_z) = -dx*dx/k_cond*(-a_e); % (i, j+1, k)
A(index, index - 1) = -dx*dx/k_cond*(-a_s-a_n); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_o); % (i, j, k+1)
end

%% near interface: index_x=N_x_1 - 1 (node ac)
vector_index_z = 2:(N_z-1);
index_x = N_x_1-1;
vector_index = vector_index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_e = x_grid_liquid(index_x+1) - x_grid_liquid(index_x); % bulk end
dx_w = x_grid_liquid(index_x) - x_grid_liquid(index_x-1); % 
k_turb_P = rho*cp*eps_H(index_x);
k_turb_E = 0; % rho*cp*eps_H(index_x+1)
k_turb_W = rho*cp*eps_H(index_x-1); 
k_turb_e = 0; % 2*k_turb_P*k_turb_E/(k_turb_P+k_turb_E);
k_turb_w = 2*k_turb_P*k_turb_W/(k_turb_P+k_turb_W);
k_cond_e = 4*k_cond_s*k_cond/(3*k_cond_s+1*k_cond); % between liquid and interface
k_cond_P = k_cond;
k_cond_w = k_cond;
k_total_P = k_cond_P + k_turb_P; % 
k_total_e = k_cond_e + k_turb_e; 
k_total_w = k_cond_w + k_turb_w; %
v_P = vz_profile(index_x);

a_w = 2*(k_total_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_total_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_total_P/dz/dz; % north, i+1
a_s = k_total_P/dz/dz + rho*cp*v_P/dz; % south, i-1
a_i = k_total_P/dy/dy; % inward, k-1
a_o = k_total_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

%% interface: index_x = N_x_1, index_z = N_z (node d)
index_x = N_x_1;
index_z = N_z;
index = index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_w = dx; % near-wall cell liquid
dx_e = dx; % near-wall cell solid
k_P = 2*k_cond_s*k_cond/(k_cond_s+k_cond);
k_w = 1/(0.5/k_cond + 0.5/k_P);
k_e = 1/(0.5/k_cond_s + 0.5/k_P);
eps_H_P = 0; % wall
eps_H_w = 0; % harmonic averaging with zero
eps_H_e = 0; % solid

a_w = 2*(k_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(index, index) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(index, index - N_z) = -dx*dx/k_cond*(-a_w); % (i, j-1, k)
A(index, index + N_z) = -dx*dx/k_cond*(-a_e); % (i, j+1, k)
A(index, index - 1) = -dx*dx/k_cond*(-a_n - a_s); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_o); % (i, j, k+1)
end

% b(index) = -(1)*dx*q_source(index_z, index_y)/k_cond;
% I_location_qbcs(index, index) = 1; % (i, j, k)

%% interface: index_x = N_x_1 (node c)
vector_index_z = 2:(N_z-1);
index_x = N_x_1;
vector_index = vector_index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_w = dx; % near-wall cell liquid
dx_e = dx; % near-wall cell solid
k_P = 2*k_cond_s*k_cond/(k_cond_s+k_cond);
k_w = 1/(0.5/k_cond + 0.5/k_P);
k_e = 1/(0.5/k_cond_s + 0.5/k_P);
eps_H_P = 0; % wall
eps_H_w = 0; % harmonic averaging with zero
eps_H_e = 0; % solid
a_w = 2*(k_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

% b(vector_index) = -(1)*dx*q_source(vector_index_z + (index_y-1)*N_z)/k_cond;
% I_location_qbcs(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = 1; % (i, j, k)

%% middle cells: index_z = N_z (node b)
vector_index_x = 2:(N_x_1-2);
index_z = N_z;
vector_index = index_z + (vector_index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

column_dx_w = x_grid_liquid(vector_index_x)' - x_grid_liquid(vector_index_x-1)'; % 
column_dx_e = x_grid_liquid(vector_index_x+1)' - x_grid_liquid(vector_index_x)'; %
column_k_turb_P = rho*cp*eps_H(vector_index_x);
column_k_turb_E = rho*cp*eps_H(vector_index_x+1);
column_k_turb_W = rho*cp*eps_H(vector_index_x-1);
column_k_turb_e = 2*column_k_turb_P.*column_k_turb_E./(column_k_turb_P+column_k_turb_E);
column_k_turb_w = 2*column_k_turb_P.*column_k_turb_W./(column_k_turb_P+column_k_turb_W);
column_k_total_P = k_cond + column_k_turb_P; % constant k_cond
column_k_total_e = k_cond + column_k_turb_e; 
column_k_total_w = k_cond + column_k_turb_w; 
column_v_P = vz_profile(vector_index_x);
column_a_w = 2*(column_k_total_w.*column_dx_e)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % west, j-1 
column_a_e = 2*(column_k_total_e.*column_dx_w)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % east, j+1
column_a_n = column_k_total_P/dz/dz; % north, i+1
column_a_s = column_k_total_P/dz/dz + rho*cp*column_v_P/dz; % south, i-1
column_a_i = column_k_total_P/dy/dy; % inward, k-1
column_a_o = column_k_total_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(column_a_w + column_a_e + column_a_n + column_a_s + column_a_i + column_a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_s-column_a_n); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_i-column_a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_i-column_a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_o);% (i, j, k+1)
end

%% middle cells (node a) 
vector_index_z = (2:(N_z-1))'; % column vector
vector_index_x = 2:(N_x_1-2); % row vector
matrix_index = repmat(vector_index_z, 1, length(vector_index_x)) + repmat(vector_index_x-1, length(vector_index_z), 1)*N_z;
vector_index = matrix_index(:) + (index_y - 1)*N_z*N_x_tot; % vectorize

matrix_dx_w = repmat(x_grid_total(vector_index_x) - x_grid_total(vector_index_x-1), length(vector_index_z), 1); % make row of dx_w and repeat in row direction
matrix_dx_e = repmat(x_grid_total(vector_index_x+1) - x_grid_total(vector_index_x), length(vector_index_z), 1);
column_dx_w = matrix_dx_w(:);
column_dx_e = matrix_dx_e(:);

matrix_k_turb_P = repmat(rho*cp*eps_H(vector_index_x)', length(vector_index_z), 1);
matrix_k_turb_E = repmat(rho*cp*eps_H(vector_index_x+1)', length(vector_index_z), 1);
matrix_k_turb_W = repmat(rho*cp*eps_H(vector_index_x-1)', length(vector_index_z), 1);
matrix_k_turb_e = 2*matrix_k_turb_P.*matrix_k_turb_E./(matrix_k_turb_P+matrix_k_turb_E);
matrix_k_turb_w = 2*matrix_k_turb_P.*matrix_k_turb_W./(matrix_k_turb_P+matrix_k_turb_W);
column_k_turb_P = matrix_k_turb_P(:);
column_k_turb_e = matrix_k_turb_e(:);
column_k_turb_w = matrix_k_turb_w(:);
column_k_total_P = k_cond + column_k_turb_P; % constant k_cond
column_k_total_e = k_cond + column_k_turb_e; 
column_k_total_w = k_cond + column_k_turb_w; 
matrix_v_P = repmat(vz_profile(vector_index_x)', length(vector_index_z), 1);
column_v_P = matrix_v_P(:);

column_a_w = 2*(column_k_total_w.*column_dx_e)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % west, j-1 
column_a_e = 2*(column_k_total_e.*column_dx_w)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % east, j+1
column_a_n = column_k_total_P/dz/dz; % north, i+1
column_a_s = column_k_total_P/dz/dz + rho*cp*column_v_P/dz; % south, i-1
column_a_i = column_k_total_P/dy/dy; % inward, k-1
column_a_o = column_k_total_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(column_a_w + column_a_e + column_a_n + column_a_s + column_a_i + column_a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_i-column_a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_i-column_a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-column_a_o);% (i, j, k+1)
end

%% assign values on matrix. note that (-1)*(1/k_liquid)*dx^2 is multiplied in equation to avoid division by small number.
% solid
%% inlet: index_z = 1
index_z = 1;
vector_index_x = N_x_1+1:N_x_tot;
vector_index = index_z + (vector_index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = 1; % (i, j, k)
% b(vector_index) = T_inlet_developed(vector_index_x);
% I_location_Tbcs(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = 1; % (i, j, k)

%% outermost wall: index_x = N_x_tot, index_z = N_z (node d')
index_z = N_z;
index_x = N_x_tot;
index = index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_w = x_grid_solid(end) - x_grid_solid(end-1); % 
dx_e = dx_w; % symmetry
k_w = k_cond_s;
k_e = k_w; % symmetry
k_P = k_cond_s;
a_w = 2*(k_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(index, index) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(index, index - N_z) = -dx*dx/k_cond*(-a_w-a_e); % (i, j-1, k)
% (i, j+1, k)
A(index, index - 1) = -dx*dx/k_cond*(-a_n - a_s); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_o); % (i, j, k+1)
end

%% outermost wall: index_x = N_x_tot (node c')
vector_index_z = 2:(N_z-1);
index_x = N_x_tot;
vector_index = vector_index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_w = x_grid_solid(end) - x_grid_solid(end-1); % 
dx_e = dx_w; % symmetry
k_w = k_cond_s;
k_e = k_w; % symmetry
k_P = k_cond_s;
a_w = 2*(k_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_w-a_e); % (i, j-1, k)
% (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

%% near interface: index_x=N_x_1 + 1, index_z = N_z (node b'd) first solid grid
index_z = N_z;
index_x = N_x_1+1;
index = index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_w = dx; % W is interface
dx_e = x_grid_solid(2) - x_grid_solid(1);
k_w = 1/(0.5/k_cond_int + 0.5/k_cond_s);
k_e = k_cond_s;
k_P = k_cond_s;
a_w = 2*(k_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(index, index) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(index, index - N_z) = -dx*dx/k_cond*(-a_w); % (i, j-1, k)
A(index, index + N_z) = -dx*dx/k_cond*(-a_e); % (i, j+1, k)
A(index, index - 1) = -dx*dx/k_cond*(-a_n - a_s); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(index, index - N_z*N_x_tot) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(index, index + N_z*N_x_tot) = -dx*dx/k_cond*(-a_o); % (i, j, k+1)
end

%% near interface: index_x=N_x_1 + 1 (node a'c)
vector_index_z = 2:(N_z-1);
index_x = N_x_1+1;
vector_index = vector_index_z + (index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

dx_w = dx; % W is interface
dx_e = x_grid_solid(2) - x_grid_solid(1);
k_w = 1/(0.5/k_cond_int + 0.5/k_cond_s);
k_e = k_cond_s;
k_P = k_cond_s;
a_w = 2*(k_w*dx_e)/(dx_e*dx_w*(dx_e+dx_w)); % west, j-1 
a_e = 2*(k_e*dx_w)/(dx_e*dx_w*(dx_e+dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(a_w + a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

%% middle cells: index_z = N_z (node b')
vector_index_x = (N_x_1+2):(N_x_tot-1);
index_z = N_z;
vector_index = index_z + (vector_index_x - 1)*N_z + (index_y - 1)*N_z*N_x_tot;

column_dx_w = x_grid_total(vector_index_x)' - x_grid_total(vector_index_x-1)';
column_dx_e = x_grid_total(vector_index_x+1)' - x_grid_total(vector_index_x)';
k_w = k_cond_s;
k_e = k_cond_s;
k_P = k_cond_s;
column_a_w = 2*(k_w.*column_dx_e)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % west, j-1 
column_a_e = 2*(k_e.*column_dx_w)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(column_a_w + column_a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s-a_n); % (i-1, j, k)
% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

%% middle cells (node a')
vector_index_z = (2:(N_z-1))'; % column vector
vector_index_x = (N_x_1+2):(N_x_tot-1); % row vector
matrix_index = repmat(vector_index_z, 1, length(vector_index_x)) + repmat(vector_index_x-1, length(vector_index_z), 1)*N_z;
vector_index = matrix_index(:) + (index_y - 1)*N_z*N_x_tot; % vectorize

matrix_dx_w = repmat(x_grid_total(vector_index_x) - x_grid_total(vector_index_x-1), length(vector_index_z), 1); % make row of dx_w and repeat in row direction
matrix_dx_e = repmat(x_grid_total(vector_index_x+1) - x_grid_total(vector_index_x), length(vector_index_z), 1);
column_dx_w = matrix_dx_w(:);
column_dx_e = matrix_dx_e(:);
k_w = k_cond_s;
k_e = k_cond_s;
k_P = k_cond_s;
column_a_w = 2*(k_w.*column_dx_e)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % west, j-1 
column_a_e = 2*(k_e.*column_dx_w)./(column_dx_e.*column_dx_w.*(column_dx_e+column_dx_w)); % east, j+1
a_n = k_P/dz/dz; % north, i+1
a_s = k_P/dz/dz; % south, i-1
a_i = k_P/dy/dy; % inward, k-1
a_o = k_P/dy/dy; % outward, k+1

A(vector_index +(vector_index-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(column_a_w + column_a_e + a_n + a_s + a_i + a_o); % (i, j, k)
A(vector_index +(vector_index-N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_w); % (i, j-1, k)
A(vector_index +(vector_index+N_z-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond.*(-column_a_e); % (i, j+1, k)
A(vector_index +(vector_index-1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_s); % (i-1, j, k)
A(vector_index +(vector_index+1-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_n);% (i+1, j, k)
if index_y == 1 % symmetry
    % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k+1)
elseif index_y == N_y % outer wall
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i-a_o); % (i, j, k-1)
    % (i, j, k+1)
else
    A(vector_index +(vector_index-N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_i); % (i, j, k-1)
    A(vector_index +(vector_index+N_z*N_x_tot-1)*N_z*N_x_tot*N_y) = -dx*dx/k_cond*(-a_o);% (i, j, k+1)
end

end % end of N_y

end