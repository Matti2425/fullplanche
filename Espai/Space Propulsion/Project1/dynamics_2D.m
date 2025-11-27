function dydt = dynamics_2D(~, y, F, c_e, mu, mode)
    r = y(1);
    vr = y(3);
    vt = y(4);
    m = y(5);
    
    % Aceleraciones gravitatorias
    ar_grav = -mu/r^2 + vt^2/r;
    at_grav = -vr*vt/r;
    
    % Aceleración del empuje
    v_mod = sqrt(vr^2 + vt^2);
    if v_mod > 1e-3
        ar_thrust = (F/m) * (vr/v_mod);
        at_thrust = (F/m) * (vt/v_mod);
    else
        ar_thrust = 0;
        at_thrust = 0;
    end
    
    dydt = [vr;
            vt/r;
            ar_grav + ar_thrust;
            at_grav + at_thrust;
            -mode * F/c_e];
end