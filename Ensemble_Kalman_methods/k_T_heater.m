% output: thermal conductivity [W/m/K] heater material
% input: temperature [K]
function [k] = k_T_heater(T)
    
    k = 0.014748*(T-273.15) + 11.1839;
    
end