function [time, H, isensors, data]=F_sensor_RSNC(U, p, maxiteration, nr)

    [n,~]=size(U);
    tic;
    % [sensors, data]=F_sensor_RSN_sub(U, p, maxiteration, nr);
    [zhat, ~, ~, ~, data] = F_sensor_RSNC_approxnt(U, p, maxiteration, nr);
    isensors = find(zhat);
    time=toc;
    [H]=F_calc_sensormatrix(p, n, isensors);
    
end
