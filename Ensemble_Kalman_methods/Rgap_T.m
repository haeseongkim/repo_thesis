% output: Rgap as temperature dependent properties (a test function)
% input: temperature [K] % output: gap thermal resistance [m^2-K/W]
function [R_gap] = Rgap_T(T_K)
    
    k_air = 0.0241 + 7.16e-5 * T_K - 1.27e-8 * T_K.^2; % thermal conductivity air, empirical correlation
    t_air = 100e-6; % [m] thickness of air gap (hyperparameter)

    e_ss = 0.8; % [-] emissivity stainless steel (hyperparameter)
    e_zro2 = 0.9; % [-] emissivity zirconia (hyperparameter)
    sigma_sb = 5.670374419e-8; % Stefan-Boltzmann constant

    h_gap_air = k_air/t_air;
    h_gap_radiation = 4*sigma_sb*T_K.^3/(1/e_ss + 1/e_zro2 - 1);
    h_gap = h_gap_air + h_gap_radiation; % [W/m^2/K]
    
    R_gap = 1./h_gap; % [m^2-K/W]

end