% velocity profile for turbulent flow : Todreas textbook eq (9.83)
% output: axial velocity
% input: radius from center [m], pipe radius [m], friction velocity [m/s], kinematic viscosity [m2/s]
function [vz_turb] = velocity_profile_turb(r, R, v_friction, nu)
    
    y_plus = (R-r)*v_friction/nu; % dimensionless distance from wall
    % follow the profile used at Kim (1989) Appendix A
    if r == R % at wall
        vz_turb = 0;
    elseif y_plus < 5% laminar sublayer
        vz_turb = v_friction*y_plus;
    elseif y_plus < 30 % buffer layer
        vz_turb = v_friction*( -3.05 + 5*log(y_plus));
    else % turbulent sublayer
%         vz_turb = v_friction*( 5.5 + 2.5*log(y_plus*1.5*(1 + r/R)/(1 + 2*(r/R)^2)) ); % Todreas
        vz_turb = v_friction*( 5.5 + 2.5*log(y_plus)); % Martinelli
    end
    
end