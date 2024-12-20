% to be used for 'display_QG_determinant_zsensors_large.m'
% find one additional sensor location which gives maximum det (smallest ellipse noise area)
function [new_indices, min_log_area_ellipse_new_sensors] = ...
    suboptimum_sensor_location_f(old_indices, T_sigma, scale_ellipsoid, H_QQ, H_GG, H_QG)

N_z = 101;
N_old = length(old_indices);
fprintf("%d\n", N_old+1)

det_new_sensors = zeros(N_z-1,1);
for id_new=1:N_z-1
    current_indices = unique([old_indices, id_new]);
    if length(current_indices) == N_old % skip iteration
    else
        det_new_sensors(id_new) = ...
            sum(H_QQ(current_indices))*sum(H_GG(current_indices)) - sum(H_QG(current_indices))^2;
    end
end

[max_det, max_idx] = max(det_new_sensors); % Find the indices of the maximum determinant
new_indices = [old_indices, max_idx];
min_log_area_ellipse_new_sensors = log10(pi*T_sigma*sqrt(scale_ellipsoid)/sqrt(max_det));

end
