function hoh = hohmannTransfer(mu, r1, r2)
    a_t = (r1 + r2)/2;
    hoh.v1 = sqrt(mu / r1); 
    hoh.v2 = sqrt(mu / r2);
    v_tp = sqrt(mu * (2/r1 - 1/a_t)); 
    v_ta = sqrt(mu * (2/r2 - 1/a_t));
    hoh.dV1 = v_tp - hoh.v1; 
    hoh.dV2 = hoh.v2 - v_ta;
    hoh.t_to_flight = pi * sqrt(a_t^3 / mu);
    hoh.a_t = a_t;
end

