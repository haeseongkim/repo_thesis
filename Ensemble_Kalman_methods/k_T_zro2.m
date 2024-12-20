% output: thermal conductivity [W/m/K] zirconia
% input: temperature [K]
function [k] = k_T_zro2(T)
    
    k = 4.8427*10^-6*(T-273.15).^2 - 0.0098*(T-273.15) + 7.5149;
    
end