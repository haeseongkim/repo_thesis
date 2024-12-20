% output: density[kg/m3]*heat capacity [J/kg/K] of zirconia
% input: temperature [K]
function [rhocp] = rhocp_T_zro2(T)
    
    rhocp = (-1.5316*10^-4*(T-273.15).^2 + 0.32*(T-273.15) + 486.0)*5700;
    
end