function xdot = ode(t, x, mu, g0, Isp, F, alpha)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    xdot = zeros(5,1);
    xdot(1) = x(3);
    xdot(2) = x(4)/x(1);
    xdot(3) = -mu/x(1)^2 + x(4)^2/x(1) + F/x(5)*sin(alpha);
    xdot(4) = -x(3)*x(4)/x(1) + F/x(5)*cos(alpha);
    xdot(5) = -F/(Isp*g0);

end