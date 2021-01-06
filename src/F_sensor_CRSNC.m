function [time, H, isensors, data]=F_sensor_CRSNC(U, p, maxiteration, nr)
% function [time_CRSN, H, sensors, NT_TOL_cal, iter]=F_sensor_CRSN(U, p, maxiteration, nr)

    [n,~]=size(U);
    tic;
    [zhat, ~, ~,~,data] = F_sensor_CRSNC_approxnt(U, p, maxiteration, nr);
    % [sensors, data]=F_sensor_CRSN_sub(U, p, maxiteration, nr);
    isensors = find(zhat);
    time=toc;
    [H]=F_calc_sensormatrix(p, n, isensors);
    
end
