function [PA_x, PA_y, PA_z] = getCoorPA(S, LA_A, LA_B)
% Gets PAs coordinates.
% S - M_B x M_A boolean support matrix - indicator if there is a sensor at each grid position.
% LA_A - locations of the sensors in the linear array A. Size M_A
% LA_B - locations of the sensors in the linear array B. Size M_B

[s_row, s_col] = find(S);

PA_x = LA_A(s_col);
PA_y = LA_B(s_row);
PA_z = zeros(size(PA_x));
