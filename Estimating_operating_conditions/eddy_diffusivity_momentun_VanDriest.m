% profile for eddy diffusivity of momentum : Van Driest
% output: eddy diffusivity of momentum [m2/s], which is used to obtain eddy diffusivity of heat
% input: radius from center [m], pipe radius [m], friction velocity [m/s], kinematic viscosity [m2/s]
function [eps_M] = eddy_diffusivity_momentun_VanDriest(r, R, v_friction, nu)
    
    A_VanDriest = 26; % 26
    K_VanDriest	= 0.4;
    
    y_plus = (R-r)*v_friction/nu; % dimensionless distance from wall
    g_y_plus = K_VanDriest*y_plus*(1-exp(-y_plus/A_VanDriest));
    eps_M = nu*0.5*(sqrt(1+4*g_y_plus*g_y_plus)-1);

end