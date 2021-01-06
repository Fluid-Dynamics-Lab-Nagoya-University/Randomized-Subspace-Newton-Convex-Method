function [ALL]=F_data_arrange1(ps, CNT, A1, A2, A3, A4)

    ALL(1:CNT, 1) = ps'        ;
    ALL(1:CNT, 2) = A1(1:CNT,1);
    ALL(1:CNT, 3) = A2(1:CNT,1);
    ALL(1:CNT, 4) = A3(1:CNT,1);
    ALL(1:CNT, 5) = A4(1:CNT,1);

end
