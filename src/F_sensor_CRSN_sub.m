function [isensors, data]=F_sensor_CRSN_sub(U, p, maxiteration, nr)
% function [isensors,data]=sensor_randomized_half_convex(U, r, p, maxiteration, nr)

    [zhat, ~, ~,~,data] = F_sensor_CRSN_approxnt(U, p, maxiteration, nr);
    % [zhat, ~, ~,~,data] = sens_sel_randomized_half_approxnt(U, p, maxiteration, nr);
    isensors = find(zhat);

end