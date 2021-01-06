function [Normalized]=F_data_normalize(ps, CNT, A1, A2, A3, A4)

    for z=1:CNT
        Normalized(z, 1) = ps(z);
        Normalized(z, 2) = A1(z,1)/A2(z,1);
        Normalized(z, 3) = A2(z,1)/A2(z,1);
        Normalized(z, 4) = A3(z,1)/A2(z,1);
        Normalized(z, 5) = A4(z,1)/A2(z,1);
    end

end


