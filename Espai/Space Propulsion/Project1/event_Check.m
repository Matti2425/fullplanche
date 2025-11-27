function [value, isterminal, direction] = event_Check(~, y, r_soi, r_planet)
    r = y(1);
    value1 = r - r_soi;
    value2 = r - (r_planet + 100e3);
    value = [value1, value2];
    isterminal = [1, 1];
    direction = [1, -1];
end