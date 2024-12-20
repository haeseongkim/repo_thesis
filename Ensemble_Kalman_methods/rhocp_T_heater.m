% output: density[kg/m3]*heat capacity [J/kg/K] of heater material
% input: temperature [K]
function [rhocp] = rhocp_T_heater(T)
    
    rhocp = (0.19963*(T-273.15) + 424.0)*8200;
    
end