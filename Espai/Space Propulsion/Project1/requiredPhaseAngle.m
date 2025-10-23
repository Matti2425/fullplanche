function phase = requiredPhaseAngle(mu, r1, r2)
    a_t = (r1 + r2)/2;
    t_trans = pi * sqrt(a_t^3 / mu);
    n2 = sqrt(mu / r2^3);
    phase = mod(pi - n2 * t_trans, 2*pi);
end