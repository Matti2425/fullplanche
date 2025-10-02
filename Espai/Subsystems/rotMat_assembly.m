function M = rotMat_assembly(axis,angle)
% axis goes 1 2 3 for x y z, respectively
switch axis
    case 1  %  X
            M = [1  0   0;
                 0  cos(angle)  -sin(angle);
                 0  sin(angle)   cos(angle)];

        case 2  % Y
            M = [ cos(angle)  0  sin(angle);
                  0  1  0;
                 -sin(angle)  0  cos(angle)];

        case 3  %  Z
            M = [ cos(angle) -sin(angle)  0;
                  sin(angle)  cos(angle)  0;
                  0           0           1];

        otherwise
            error('Axis error')

end